'use strict';

// Viewer for viral-sim: runs a simulation through the local server
// (/api/run) or loads a JSON file from `viral-sim run --json`, and draws
// what was computed, layer by layer. Studies (many runs, results as curves)
// run in the background on the server (/api/study/...) and are polled.
// Plain canvas, no libraries.

const $ = (id) => document.getElementById(id);
const form = $('controls');
const F = (name) => form.elements[name];
const studyForm = $('studyForm');
const SF = (name) => studyForm.elements[name];

let data = null;				// the last run
let layer = 0;					// layer shown (0 = initial state)
let playing = false;
let setsMode = 'layer';			// split shown in the problem picture: layer, final, local
let skView = 'j';				// SK heatmap: j or frustration
let mode = 'run';				// run or studies
let study = null;				// the study shown (partial or final)
let studyRunning = false;
let sweepView = 'grid';			// sweep heatmap: grid, grid_local, grid_energy
const view = { yaw: -0.65, pitch: 0.38 };	// sphere orientation
const hover = {};				// canvas id -> hovered index

const DEFAULTS = {
	partition: { n: 24, layers: 100, dt: 0.5 },
	sk: { n: 300, layers: 200, dt: 0.3 },
	'sk-gauss': { n: 300, layers: 200, dt: 0.3 },
	maxcut: { n: 100, layers: 200, dt: 0.3 },
	chain: { n: 100, layers: 1000, dt: 0.3 },
	hopfield: { n: 400, layers: 100, dt: 0.3 },
};

const MAX_N = { partition: 2000000, sk: 12000, 'sk-gauss': 12000, maxcut: 12000, chain: 1000000, hopfield: 200000 };

// extra inputs, and the problem that uses each
const EXTRAS = { degree: 'maxcut', alpha: 'chain', patterns: 'hopfield', noise: 'hopfield' };

const STUDY_PRESETS = {
	scaling: {
		quick: { problem: 'sk', sizes: '250,500,1000,2000,4000,8000', seeds: 10, layers: 200, dt: 0.3 },
		big: { problem: 'sk', sizes: '500,1000,2000,4000,8000,12000', seeds: 20, layers: 500, dt: 0.5 },
	},
	sweep: {
		quick: { problem: 'sk', n: 300, layers_list: '25,50,100,200,400,800', dt_list: '0.1,0.2,0.3,0.5,0.8,1.2', seeds: 5 },
		big: { problem: 'maxcut', n: 2000, degree: 3, layers_list: '25,50,100,200,400,800,1600', dt_list: '0.1,0.2,0.3,0.4,0.5,0.7,0.9,1.2', seeds: 5 },
	},
	squeezing: {
		quick: { n: 1000, trajectories: 2000, times: 60, t_max: '' },
		big: { n: 100000, trajectories: 20000, times: 60, t_max: '' },
	},
};


// ---------------------------------------------------------------
// colors, numbers, canvas
// ---------------------------------------------------------------

function readTheme() {
	const s = getComputedStyle(document.documentElement);
	const g = (n) => s.getPropertyValue(n).trim();
	return {
		text: g('--text'), text2: g('--text-2'), muted: g('--muted'), grid: g('--grid'),
		axis: g('--axis'), surface: g('--surface'), a: g('--set-a'), b: g('--set-b'),
		c3: g('--series-3'), c7: g('--series-7'), pos: g('--pos'), neg: g('--neg'), mid: g('--mid'),
	};
}
let T = readTheme();

function rgb(hex) {
	let h = hex.replace('#', '');
	if (h.length === 3) h = h.split('').map((c) => c + c).join('');
	const n = parseInt(h, 16);
	return [(n >> 16) & 255, (n >> 8) & 255, n & 255];
}

function mixRgb(c1, c2, t) {
	const a = rgb(c1), b = rgb(c2);
	return a.map((v, i) => Math.round(v + (b[i] - v) * t));
}

const css = (c) => `rgb(${c[0]},${c[1]},${c[2]})`;

function fmt(v, digits = 4) {
	if (v === null || v === undefined || !isFinite(v)) return '—';
	if (v === 0) return '0';
	const a = Math.abs(v);
	if (a >= 1e6 || a < 1e-3) {
		const e = Math.floor(Math.log10(a));
		return `${Number((v / Math.pow(10, e)).toPrecision(3))}×${pow10(e)}`;
	}
	return String(Number(v.toPrecision(digits)));
}

const pct = (v, d = 2) => (v === null || v === undefined || !isFinite(v) ? '—' : (100 * v).toFixed(d) + '%');

function fmtTime(s) {
	if (s === null || s === undefined || !isFinite(s)) return '—';
	if (s < 1e-6) return (s * 1e9).toFixed(0) + ' ns';
	if (s < 1e-3) return (s * 1e6).toFixed(1) + ' µs';
	if (s < 1) return (s * 1e3).toFixed(2) + ' ms';
	if (s < 120) return s.toFixed(s < 10 ? 2 : 1) + ' s';
	return `${Math.floor(s / 60)} min ${Math.round(s % 60)} s`;
}

const SUP = { '-': '⁻', 0: '⁰', 1: '¹', 2: '²', 3: '³', 4: '⁴', 5: '⁵', 6: '⁶', 7: '⁷', 8: '⁸', 9: '⁹' };
const pow10 = (k) => (k === 0 ? '1' : k === 1 ? '10' : '10' + String(k).split('').map((c) => SUP[c]).join(''));

function niceTicks(lo, hi, count) {
	if (!(hi > lo)) return [lo];
	const raw = (hi - lo) / Math.max(1, count);
	const mag = Math.pow(10, Math.floor(Math.log10(raw)));
	const e = raw / mag;
	const step = (e >= 7.5 ? 10 : e >= 3.5 ? 5 : e >= 1.5 ? 2 : 1) * mag;
	const out = [];
	for (let t = Math.ceil(lo / step) * step; t <= hi + step * 1e-9; t += step)
		out.push(Math.abs(t) < step * 1e-9 ? 0 : Number(t.toPrecision(12)));
	return out;
}

// Sizes the canvas for the screen's pixel ratio and clears it.
function setup(cv) {
	const dpr = window.devicePixelRatio || 1;
	const w = cv.clientWidth, h = cv.clientHeight;
	const W = Math.max(1, Math.round(w * dpr)), H = Math.max(1, Math.round(h * dpr));
	if (cv.width !== W || cv.height !== H) { cv.width = W; cv.height = H; }
	const ctx = cv.getContext('2d');
	ctx.setTransform(dpr, 0, 0, dpr, 0, 0);
	ctx.clearRect(0, 0, w, h);
	ctx.font = '11px ' + getComputedStyle(document.body).fontFamily;
	return { ctx, w, h };
}

function line(ctx, x1, y1, x2, y2) {
	ctx.beginPath();
	ctx.moveTo(x1, y1);
	ctx.lineTo(x2, y2);
	ctx.stroke();
}

function dot(ctx, x, y, r, color) {
	ctx.beginPath();
	ctx.arc(x, y, r, 0, 2 * Math.PI);
	ctx.fillStyle = color;
	ctx.fill();
	ctx.lineWidth = 2;
	ctx.strokeStyle = T.surface;
	ctx.stroke();
}

// text with a surface-colored halo, so it stays readable over lines
function haloText(ctx, text, x, y, color) {
	ctx.lineWidth = 3;
	ctx.strokeStyle = T.surface;
	ctx.lineJoin = 'round';
	ctx.strokeText(text, x, y);
	ctx.fillStyle = color;
	ctx.fillText(text, x, y);
}

function key(color, label, kind = '') {
	return `<span class="key"><span class="sw ${kind}" style="${kind === 'dash' ? '' : 'background:' + color}"></span>${label}</span>`;
}

// a 2D grid of values in [-1, 1] as an image: -1 light, +1 dark
function drawPixels(ctx, vals, rows, cols, x, y, size) {
	const img = new ImageData(cols, rows);
	const lo = rgb(T.mid), hi = rgb(T.text);
	for (let k = 0; k < rows * cols; k++) {
		const o = 4 * k;
		if (k >= vals.length) { img.data[o + 3] = 0; continue; }
		const t = (Math.max(-1, Math.min(1, vals[k])) + 1) / 2;
		for (let c = 0; c < 3; c++) img.data[o + c] = Math.round(lo[c] + (hi[c] - lo[c]) * t);
		img.data[o + 3] = 255;
	}
	const off = document.createElement('canvas');
	off.width = cols; off.height = rows;
	off.getContext('2d').putImageData(img, 0, 0);
	ctx.imageSmoothingEnabled = false;
	const cell = size / Math.max(rows, cols);
	ctx.drawImage(off, x, y, cols * cell, rows * cell);
	ctx.strokeStyle = T.axis;
	ctx.lineWidth = 1;
	ctx.strokeRect(x - 0.5, y - 0.5, cols * cell + 1, rows * cell + 1);
}


// ---------------------------------------------------------------
// tooltip and hover
// ---------------------------------------------------------------

const tip = $('tooltip');

function showTip(html, x, y) {
	tip.innerHTML = html;
	tip.hidden = false;
	const r = tip.getBoundingClientRect();
	let left = x + 14, top = y + 14;
	if (left + r.width > window.innerWidth - 8) left = x - r.width - 14;
	if (top + r.height > window.innerHeight - 8) top = y - r.height - 14;
	tip.style.left = Math.max(4, left) + 'px';
	tip.style.top = Math.max(4, top) + 'px';
}

function hideTip() { tip.hidden = true; }

// Charts store canvas._geom = { pick(px, py), tip(i), render(), click(i) }.
function attachHover(cv) {
	cv.addEventListener('mousemove', (e) => {
		const g = cv._geom;
		if (!g || drag) return;
		const r = cv.getBoundingClientRect();
		const i = g.pick(e.clientX - r.left, e.clientY - r.top);
		if (i !== hover[cv.id]) { hover[cv.id] = i; g.render(); }
		if (i === null || i === undefined) hideTip();
		else showTip(cv._geom.tip(i), e.clientX, e.clientY);
	});
	cv.addEventListener('mouseleave', () => {
		hover[cv.id] = null;
		hideTip();
		if (cv._geom) cv._geom.render();
	});
	cv.addEventListener('click', () => {
		const g = cv._geom, i = hover[cv.id];
		if (g && g.click && i !== null && i !== undefined) g.click(i);
	});
}

// dashed reference lines; ones closer than 4 px share a line and a label
function drawRefLines(ctx, refs, Y, x0, x1, top) {
	const groups = [];
	for (const l of refs) {
		const y = Y(l.value);
		if (!isFinite(y)) continue;
		const g = groups.find((g) => Math.abs(g.y - y) < 4);
		if (g) g.names.push(l.label); else groups.push({ y, names: [l.label] });
	}
	ctx.save();
	ctx.setLineDash([4, 4]);
	ctx.strokeStyle = T.muted;
	ctx.lineWidth = 1;
	for (const g of groups) line(ctx, x0, g.y, x1, g.y);
	ctx.restore();
	return () => {
		ctx.font = '11px ' + getComputedStyle(document.body).fontFamily;
		ctx.textAlign = 'left';
		ctx.textBaseline = 'bottom';
		const labels = groups.map((g) => ({
			y: g.y,
			text: g.names.length <= 2 ? g.names.join(' / ') : `${g.names[0]} and ${g.names.length - 1} others (nearly equal)`,
		})).sort((a, b) => b.y - a.y);
		let last = Infinity;
		for (const lb of labels) {
			const y = Math.min(lb.y - 2, last - 12);
			last = y;
			if (y < top + 10) continue;
			haloText(ctx, lb.text, x0 + 6, y, T.text2);
		}
	};
}


// ---------------------------------------------------------------
// line chart: series over layers, dashed baselines, layer marker
// ---------------------------------------------------------------

// o: { series: [{name, values, color, x0}], hlines: [{label, value}],
//      logY, floor, marker, render, click, xName }
function lineChart(cv, o) {
	const { ctx, w, h } = setup(cv);
	const pad = { l: 58, r: 14, t: 10, b: 26 };
	const pw = w - pad.l - pad.r, ph = h - pad.t - pad.b;
	const xmin = o.xmin ?? 0;
	const xmax = Math.max(xmin + 1, ...o.series.map((s) => (s.x0 || 0) + s.values.length - 1));
	const floor = o.floor ?? 1e-12;
	const ty = o.logY ? (v) => Math.log10(Math.max(v, floor)) : (v) => v;

	let lo = Infinity, hi = -Infinity;
	const take = (v) => {
		if (v === null || v === undefined || !isFinite(v)) return;
		const t = ty(v);
		if (t < lo) lo = t;
		if (t > hi) hi = t;
	};
	o.series.forEach((s) => s.values.forEach(take));
	(o.hlines || []).forEach((l) => { if (!(o.logY && !(l.value > 0))) take(l.value); });
	if (!isFinite(lo)) { lo = 0; hi = 1; }
	if (hi - lo < 1e-12) { lo -= 0.5; hi += 0.5; }
	const m = (hi - lo) * 0.06;
	lo -= m; hi += m;

	const X = (i) => pad.l + (i - xmin) / (xmax - xmin) * pw;
	const Yt = (t) => pad.t + (hi - t) / (hi - lo) * ph;
	const Y = (v) => Yt(ty(v));

	// grid and y labels
	ctx.textAlign = 'right';
	ctx.textBaseline = 'middle';
	ctx.lineWidth = 1;
	const yticks = o.logY
		? Array.from({ length: Math.floor(hi) - Math.ceil(lo) + 1 }, (_, k) => Math.ceil(lo) + k)
		: niceTicks(lo, hi, 5);
	const step = o.logY ? Math.max(1, Math.ceil(yticks.length / 6)) : 1;
	yticks.forEach((t, k) => {
		if (k % step) return;
		const y = Yt(t);
		ctx.strokeStyle = T.grid;
		line(ctx, pad.l, y, pad.l + pw, y);
		ctx.fillStyle = T.muted;
		ctx.fillText(o.logY ? pow10(t) : fmt(t, 3), pad.l - 6, y);
	});

	// x labels
	ctx.textAlign = 'center';
	ctx.textBaseline = 'top';
	ctx.fillStyle = T.muted;
	for (const t of niceTicks(xmin, xmax, Math.max(2, Math.min(8, Math.floor(pw / 60))))) {
		if (t !== Math.round(t)) continue;
		ctx.fillText(String(t), X(t), pad.t + ph + 6);
	}
	ctx.strokeStyle = T.axis;
	line(ctx, pad.l, pad.t + ph + 0.5, pad.l + pw, pad.t + ph + 0.5);

	const refs = (o.hlines || []).filter((l) => l.value !== null && l.value !== undefined && !(o.logY && !(l.value > 0)));
	const labelRefs = drawRefLines(ctx, refs.map((l) => ({ label: l.label, value: l.value })), Y, pad.l, pad.l + pw, pad.t);

	// series, the ones marked `under` first
	for (const s of [...o.series.filter((s) => s.under), ...o.series.filter((s) => !s.under)]) {
		ctx.strokeStyle = s.color;
		ctx.lineWidth = 2;
		ctx.lineJoin = 'round';
		ctx.beginPath();
		let on = false;
		s.values.forEach((v, k) => {
			if (v === null || v === undefined || !isFinite(v)) { on = false; return; }
			const x = X((s.x0 || 0) + k), y = Y(v);
			if (on) ctx.lineTo(x, y); else { ctx.moveTo(x, y); on = true; }
		});
		ctx.stroke();
	}
	labelRefs();

	// current layer and hover
	const marks = (i, color) => {
		const x = X(i);
		ctx.strokeStyle = color;
		ctx.lineWidth = 1;
		line(ctx, x, pad.t, x, pad.t + ph);
		for (const s of o.series) {
			const v = s.values[i - (s.x0 || 0)];
			if (v !== null && v !== undefined && isFinite(v)) dot(ctx, x, Y(v), 4, s.color);
		}
	};
	if (o.marker !== null && o.marker !== undefined && o.marker >= xmin && o.marker <= xmax) marks(o.marker, T.text2);
	const hv = hover[cv.id];
	if (hv !== null && hv !== undefined) marks(hv, T.muted);

	cv._geom = {
		render: o.render,
		click: o.click,
		pick: (px) => {
			if (px < pad.l - 10 || px > pad.l + pw + 10) return null;
			return Math.max(xmin, Math.min(xmax, Math.round(xmin + (px - pad.l) / pw * (xmax - xmin))));
		},
		tip: (i) => `<b>${o.xName || 'Layer'} ${i}</b>` + o.series.map((s) => {
			const v = s.values[i - (s.x0 || 0)];
			return v === undefined ? '' : `<div><span class="sw line" style="background:${s.color}"></span>${s.name}: ${fmt(v, 6)}</div>`;
		}).join('') + (o.click ? '<div class="note">click to show this layer</div>' : ''),
	};
}


// ---------------------------------------------------------------
// x-y chart for the studies: points, error bars, fitted lines
// ---------------------------------------------------------------

// o: { series: [{name, color, xs, ys, err, line, dots, r, alpha}], hlines: [{label, value}],
//      fits: [{color, a, b, x0, x1, label}], marks: [{x, y, color, label}],
//      xmin, xmax, yFmt, xFmt, topTicks: [{x, label}], xTitle, pickXs, tip(i), render }
function xyChart(cv, o) {
	const { ctx, w, h } = setup(cv);
	const pad = { l: 58, r: 16, t: o.topTicks ? 30 : 10, b: o.xTitle ? 42 : 26 };
	const pw = w - pad.l - pad.r, ph = h - pad.t - pad.b;
	const xsAll = o.series.flatMap((s) => s.xs);
	let xmin = o.xmin ?? Math.min(...xsAll), xmax = o.xmax ?? Math.max(...xsAll);
	if (!isFinite(xmin) || !isFinite(xmax)) { xmin = 0; xmax = 1; }
	if (xmax - xmin < 1e-15) { xmin -= 0.5; xmax += 0.5; }

	let lo = Infinity, hi = -Infinity;
	const take = (v) => { if (Number.isFinite(v)) { lo = Math.min(lo, v); hi = Math.max(hi, v); } };
	for (const s of o.series) s.ys.forEach((y, k) => {
		const e = s.err ? s.err[k] || 0 : 0;
		take(y - e); take(y + e);
	});
	(o.hlines || []).forEach((l) => take(l.value));
	(o.fits || []).forEach((f) => { take(f.a + f.b * f.x0); take(f.a + f.b * f.x1); });
	if (!isFinite(lo)) { lo = 0; hi = 1; }
	if (hi - lo < 1e-12) { lo -= 0.5; hi += 0.5; }
	const m = (hi - lo) * 0.07;
	lo -= m; hi += m;

	const X = (x) => pad.l + (x - xmin) / (xmax - xmin) * pw;
	const Y = (y) => pad.t + (hi - y) / (hi - lo) * ph;
	const yFmt = o.yFmt || ((v) => fmt(v, 3));
	const xFmt = o.xFmt || ((v) => fmt(v, 3));

	ctx.textAlign = 'right';
	ctx.textBaseline = 'middle';
	ctx.lineWidth = 1;
	for (const t of niceTicks(lo, hi, 5)) {
		const y = Y(t);
		ctx.strokeStyle = T.grid;
		line(ctx, pad.l, y, pad.l + pw, y);
		ctx.fillStyle = T.muted;
		ctx.fillText(yFmt(t), pad.l - 6, y);
	}
	ctx.textAlign = 'center';
	ctx.textBaseline = 'top';
	for (const t of niceTicks(xmin, xmax, Math.max(2, Math.min(7, Math.floor(pw / 70))))) {
		ctx.fillStyle = T.muted;
		ctx.fillText(xFmt(t), X(t), pad.t + ph + 6);
	}
	ctx.strokeStyle = T.axis;
	line(ctx, pad.l, pad.t + ph + 0.5, pad.l + pw, pad.t + ph + 0.5);
	if (o.xTitle) {
		ctx.fillStyle = T.text2;
		ctx.fillText(o.xTitle, pad.l + pw / 2, pad.t + ph + 22);
	}
	if (o.topTicks) {
		ctx.textBaseline = 'bottom';
		let last = -Infinity;
		for (const tk of [...o.topTicks].sort((a, b) => X(b.x) - X(a.x))) {
			const x = X(tk.x);
			ctx.strokeStyle = T.axis;
			line(ctx, x, pad.t - 4, x, pad.t);
			const wdt = ctx.measureText(tk.label).width;
			if (Math.abs(x - last) < wdt + 8) continue;
			last = x;
			ctx.fillStyle = T.muted;
			ctx.fillText(tk.label, Math.min(pad.l + pw - wdt / 2, Math.max(pad.l + wdt / 2, x)), pad.t - 6);
		}
		ctx.strokeStyle = T.axis;
		line(ctx, pad.l, pad.t - 0.5, pad.l + pw, pad.t - 0.5);
	}

	const labelRefs = drawRefLines(ctx, (o.hlines || []).map((l) => ({ label: l.label, value: l.value })), Y, pad.l, pad.l + pw, pad.t);

	for (const f of o.fits || []) {
		ctx.save();
		ctx.setLineDash([6, 4]);
		ctx.strokeStyle = f.color;
		ctx.lineWidth = 1.5;
		line(ctx, X(f.x0), Y(f.a + f.b * f.x0), X(f.x1), Y(f.a + f.b * f.x1));
		ctx.restore();
	}

	for (const s of o.series) {
		if (s.line) {
			ctx.strokeStyle = s.color;
			ctx.lineWidth = 2;
			ctx.lineJoin = 'round';
			ctx.beginPath();
			let on = false;
			s.ys.forEach((y, k) => {
				if (!Number.isFinite(y)) { on = false; return; }
				if (on) ctx.lineTo(X(s.xs[k]), Y(y)); else { ctx.moveTo(X(s.xs[k]), Y(y)); on = true; }
			});
			ctx.stroke();
		}
		if (s.err) {
			ctx.strokeStyle = s.color;
			ctx.lineWidth = 1.5;
			s.ys.forEach((y, k) => {
				const e = s.err[k];
				if (!Number.isFinite(y) || !(e > 0)) return;
				const x = X(s.xs[k]);
				line(ctx, x, Y(y - e), x, Y(y + e));
				line(ctx, x - 4, Y(y - e), x + 4, Y(y - e));
				line(ctx, x - 4, Y(y + e), x + 4, Y(y + e));
			});
		}
		if (s.dots) {
			ctx.globalAlpha = s.alpha ?? 1;
			s.ys.forEach((y, k) => {
				if (!Number.isFinite(y)) return;
				if ((s.r ?? 4) < 3) {
					ctx.beginPath();
					ctx.arc(X(s.xs[k]), Y(y), s.r, 0, 2 * Math.PI);
					ctx.fillStyle = s.color;
					ctx.fill();
				} else dot(ctx, X(s.xs[k]), Y(y), s.r ?? 4, s.color);
			});
			ctx.globalAlpha = 1;
		}
	}

	for (const f of o.fits || []) {
		if (!f.label) continue;
		const x = X(f.x0), y = Y(f.a + f.b * f.x0);
		ctx.beginPath();
		ctx.arc(x, y, 5, 0, 2 * Math.PI);
		ctx.fillStyle = T.surface;
		ctx.fill();
		ctx.strokeStyle = f.color;
		ctx.lineWidth = 2;
		ctx.stroke();
		ctx.textAlign = 'left';
		ctx.textBaseline = f.below ? 'top' : 'bottom';
		haloText(ctx, f.label, x + 8, y + (f.below ? 4 : -4), f.color);
	}
	for (const mk of o.marks || []) {
		const x = X(mk.x), y = Y(mk.y);
		if (!isFinite(y)) continue;
		ctx.beginPath();
		ctx.arc(x, y, 7, 0, 2 * Math.PI);
		ctx.strokeStyle = mk.color;
		ctx.lineWidth = 2;
		ctx.stroke();
		ctx.textAlign = x > pad.l + pw * 0.6 ? 'right' : 'left';
		ctx.textBaseline = mk.above ? 'bottom' : 'top';
		haloText(ctx, mk.label, x + (ctx.textAlign === 'left' ? 10 : -10), y + (mk.above ? -8 : 8), mk.color);
	}
	labelRefs();

	const hv = hover[cv.id];
	if (o.pickXs && hv !== null && hv !== undefined) {
		ctx.strokeStyle = T.muted;
		ctx.lineWidth = 1;
		const x = X(o.pickXs[hv]);
		line(ctx, x, pad.t, x, pad.t + ph);
	}

	cv._geom = o.pickXs ? {
		render: o.render,
		pick: (px) => {
			if (px < pad.l - 12 || px > pad.l + pw + 12) return null;
			let best = null, bd = Infinity;
			o.pickXs.forEach((x, k) => {
				const d = Math.abs(X(x) - px);
				if (d < bd) { bd = d; best = k; }
			});
			return bd < 40 ? best : null;
		},
		tip: o.tip,
	} : null;
}


// ---------------------------------------------------------------
// what is being solved
// ---------------------------------------------------------------

function renderAbout() {
	const p = data.params, kind = data.problem.kind, n = data.problem.n, pr = data.problem;
	let html = '';
	if (kind === 'partition') {
		html += `<p><b>Number partitioning.</b> Split the ${n.toLocaleString()} numbers a<sub>i</sub>${pr.weights ? ' (the bars)' : ''} into two sets whose sums are as equal as possible. Spin i's sign s<sub>i</sub> = ±1 says which set number i goes in, so the imbalance is |Σ a<sub>i</sub> s<sub>i</sub>|. Making it small is the same as minimizing the energy</p>
			<p class="formula">E = Σ<sub>i&lt;j</sub> a<sub>i</sub> a<sub>j</sub> s<sub>i</sub> s<sub>j</sub> = ((Σ a<sub>i</sub> s<sub>i</sub>)² − Σ a<sub>i</sub>²) / 2.</p>
			<p>The coupling J<sub>ij</sub> = a<sub>i</sub> a<sub>j</sub> has rank 1, which is what <code>VIRAL_n_spin</code> computes on the FPGA: every spin feels the others only through S = Σ a<sub>j</sub> z<sub>j</sub>, so a layer costs N steps, not N².</p>`;
	} else if (kind === 'sk' || kind === 'sk-gauss') {
		const how = kind === 'sk' ? '±1/√N' : 'Gaussian, with variance 1/N';
		html += `<p><b>Sherrington–Kirkpatrick spin glass.</b> Every pair of the ${n} spins is coupled with a random strength J<sub>ij</sub> (${how}), so J is a dense ${n}×${n} matrix. The goal is the configuration with the lowest energy</p>
			<p class="formula">E = Σ<sub>i&lt;j</sub> J<sub>ij</sub> s<sub>i</sub> s<sub>j</sub>.</p>
			<p>The couplings contradict each other (frustration), so no configuration satisfies them all; the best energy approaches −0.763 per spin as N grows (Parisi). A dense J doesn't fit <code>VIRAL_n_spin</code>: this is the computation of the proposed larger FPGA design, run here on the CPU.</p>`;
	} else if (kind === 'maxcut') {
		const m = data.result.maxcut;
		html += `<p><b>Max-Cut.</b> The ${n} vertices of a random graph, each with ${m.degree} neighbours (${m.edges.toLocaleString()} edges), are split into two sides. An edge is cut when its two ends are on different sides, and the goal is to cut as many edges as possible. With J = 1 on every edge the energy is</p>
			<p class="formula">E = Σ<sub>edges</sub> s<sub>i</sub> s<sub>j</sub> = edges − 2 × cut,</p>
			<p>so the lowest energy is the largest cut. This is the standard benchmark for QAOA: on large random 3-regular graphs, QAOA with one layer cuts 69.24% of the edges and with two layers 75.59%.</p>
			<p>J is sparse: each spin has only ${m.degree} partners, so a dense J wastes almost all of its memory and multiply-adds. A sparse-J design (like the p-bit machines) would need only N × ${m.degree} multiply-adds per layer. Here J is stored dense.</p>`;
	} else if (kind === 'chain') {
		html += `<p><b>Long-range Ising chain.</b> ${n.toLocaleString()} spins in a row, every pair coupled antiferromagnetically with strength J = 1/r<sup>α</sup>, where r is their distance (α = ${fmt(pr.alpha, 3)}). This is the Ising chain of trapped-ion quantum simulators (α between 0 and 3) and of dipolar molecules (α = 3). Every coupling prefers the two spins to point opposite ways, so the ground state is perfect alternation ↑↓↑↓.</p>
			<p>A slow anneal reaches it. A fast one leaves defects: stretches of the chain alternate in the opposite phase, separated by domain walls where two neighbours point the same way. Fewer layers leave more defects (the Kibble–Zurek effect); try 50 or 150 layers.</p>
			<p>J depends only on distance, so an FPGA wouldn't store it: it would look 1/r<sup>α</sup> up in a small table, N² multiply-adds per layer and no J memory. The CPU gets the same fields from an FFT in N log N steps, so for long chains the CPU has the better algorithm.</p>`;
	} else if (kind === 'hopfield') {
		const hp = data.result.hopfield;
		const pic = pr.shape ? `, drawn as ${pr.shape[0]}×${pr.shape[1]} images` : '';
		html += `<p><b>Hopfield associative memory.</b> ${hp.patterns} random patterns ξ<sup>μ</sup> of ±1 (${n} spins each${pic}) are stored in the couplings</p>
			<p class="formula">J = −(1/N) Σ<sub>μ</sub> ξ<sup>μ</sup> ξ<sup>μ</sup>ᵀ,</p>
			<p>which makes every stored pattern, and its inverse, a low-energy state. J has rank P = ${hp.patterns}: between the rank 1 that <code>VIRAL_n_spin</code> computes and a full J, it needs P sums per layer instead of one.</p>
			<p>Classic recall starts from a noisy cue (pattern 1 with ${Math.round(100 * (pr.noise ?? p.noise))}% of its spins flipped) and flips spins downhill; it works up to a load P/N ≈ 0.138. Annealing is different: it looks for the lowest energy overall, so it finds a stored pattern only below P/N ≈ 0.05, where the stored patterns are the lowest-energy states (Amit, Gutfreund and Sompolinsky). Between 0.05 and 0.138 they are only local minima, and the anneal ends in a mixture instead. This run: P/N = <b>${fmt(hp.load, 3)}</b>.</p>
			<p>The start is tilted toward the cue, but the anneal mostly ignores it: it may find another stored pattern, or an inverted one.</p>`;
	}
	const start = p.init === 'x'
		? (kind === 'hopfield' ? `along +x, with z tilted by ${p.tilt} toward the cue` : `along +x, with z tilted by up to ±${p.tilt} so the symmetry breaks`)
		: 'in random directions';
	html += `<p><b>What is computed.</b> Each spin is a classical unit vector (x, y, z), the mean-field picture of a qubit. All start ${start}. In layer k:</p>
		<ol>
			<li><b>Cost step:</b> spin i turns about z by an angle proportional to γ<sub>k</sub> h<sub>i</sub>, where h<sub>i</sub> = Σ<sub>j</sub> J<sub>ij</sub> z<sub>j</sub> is the field from the other spins.</li>
			<li><b>Mixer step:</b> every spin turns about x by 2β<sub>k</sub>.</li>
		</ol>
		<p>${p.schedule === 'ramp'
			? `Over the ${p.layers} layers γ rises and β falls (a linear ramp of size dt = ${p.dt}), so the spins move from the mixer's state toward low energy.`
			: 'This run uses the hardware\'s 4-entry angle_lut, a test pattern rather than an optimizing schedule.'}
		At the end each spin is rounded to s<sub>i</sub> = sign(z<sub>i</sub>). A local search then polishes that answer, and classical methods give a comparison.</p>`;
	const engine = {
		partition: 'Double-precision spins and weights.',
		sk: 'Double-precision spins; the dense J is stored as 8-bit integers, as an FPGA would store it.',
		'sk-gauss': 'Double-precision spins; the dense J is stored as 8-bit integers, as an FPGA would store it.',
		maxcut: 'Double-precision spins; J (1 on each edge) is stored as a dense 8-bit matrix.',
		chain: 'Double-precision spins; the fields come from an FFT of 1/r^α, exactly.',
		hopfield: 'Double-precision spins; J is kept as the P patterns, so a layer costs N × P steps.',
	}[kind];
	html += `<p><span class="tag">${p.engine === 'hw' ? 'engine: hw' : 'engine: float'}</span>${p.engine === 'hw'
		? 'Bit-exact model of the Verilog: 16-bit spins, CORDIC rotations, the same integers as <code>VIRAL_n_spin</code>. The floating-point result is drawn alongside.'
		: engine}</p>`;
	$('about').innerHTML = html;
}


// ---------------------------------------------------------------
// problem picture
// ---------------------------------------------------------------

function finalAssignment(local) {
	const r = data.result;
	if (r.partition) return (local ? r.partition.assignment_local : r.partition.assignment) || null;
	return (local ? r.assignment_local : r.assignment) || null;
}

const allShown = () => data.shown.length === data.problem.n;

// signs of z at layer k, in `shown` order
function layerSigns(k) {
	const s = data.layers[k].spins, out = new Array(data.shown.length);
	for (let m = 0; m < out.length; m++) out[m] = s[3 * m + 2] < 0 ? -1 : 1;
	return out;
}

// z of every spin at layer k (all spins shown)
function layerZ(k) {
	const s = data.layers[k].spins, out = new Array(data.shown.length);
	for (let m = 0; m < out.length; m++) out[m] = s[3 * m + 2];
	return out;
}

// signs of every spin for the split being shown, or null
function currentSigns() {
	if (setsMode === 'layer') return allShown() ? layerSigns(layer) : null;
	return finalAssignment(setsMode === 'local');
}

function updateToggles() {
	const kind = data.problem.kind;
	for (const b of $('setsMode').querySelectorAll('button')) {
		const m = b.dataset.mode;
		b.disabled = (m === 'layer' && !allShown()) || (m !== 'layer' && !finalAssignment(m === 'local'));
		b.classList.toggle('on', m === setsMode);
	}
	$('setsMode').hidden = [...$('setsMode').querySelectorAll('button')].every((b) => b.disabled);
	$('skView').hidden = !(kind === 'sk' || kind === 'sk-gauss') || !data.problem.j;
	for (const b of $('skView').querySelectorAll('button')) b.classList.toggle('on', b.dataset.view === skView);
}

// resets the parts of the problem card that only some problems use
function problemCard(title) {
	$('problemTitle').textContent = title;
	$('problemCanvas').hidden = false;
	$('problemCanvas').style.height = '';
	$('problemCanvas2').hidden = true;
	$('problemCaption2').hidden = true;
	$('problemLegend').innerHTML = '';
	$('problemNote').innerHTML = '';
	$('problemCanvas2')._geom = null;
}

function renderProblem() {
	updateToggles();
	const kind = data.problem.kind;
	if (kind === 'partition') renderPartition();
	else if (kind === 'maxcut') renderMaxCut();
	else if (kind === 'chain') renderChain();
	else if (kind === 'hopfield') renderHopfield();
	else renderSK();
}

function renderPartition() {
	const cv = $('problemCanvas'), w = data.problem.weights, part = data.result.partition;
	problemCard('The numbers and the split');
	if (!w) {
		cv.hidden = true;
		$('problemCaption').textContent = `The numbers are listed only for N ≤ 20,000, so there are no bars here. The final mean-field split of all ${data.problem.n.toLocaleString()}:`;
		$('problemNote').innerHTML = totalsHtml(part.sum_a, part.sum_b, part.residual, part);
		return;
	}
	const signs = currentSigns() || finalAssignment(false);
	const order = w.map((_, i) => i).sort((i, j) => w[j] - w[i]);
	const MAX = 400;
	const picked = order.length > MAX ? Array.from({ length: MAX }, (_, k) => order[Math.floor(k * order.length / MAX)]) : order;
	$('problemCaption').innerHTML = `Each bar is one number a<sub>i</sub>, largest first${picked.length < order.length ? ` (${picked.length} of ${order.length} shown)` : ''}; the color is the set its spin puts it in. ${
		setsMode === 'layer' ? `Split at layer ${layer}: the sign of each z now.` : setsMode === 'local' ? 'Final split after the local search.' : 'Final split: the signs of z after the last layer.'}`;
	cv.style.height = '210px';
	const { ctx, w: W, h: H } = setup(cv);
	const pad = { l: 40, r: 6, t: 8, b: 8 };
	const pw = W - pad.l - pad.r, ph = H - pad.t - pad.b;
	const top = Math.max(...niceTicks(0, w[order[0]], 4));
	ctx.textAlign = 'right';
	ctx.textBaseline = 'middle';
	ctx.lineWidth = 1;
	for (const t of niceTicks(0, top, 4)) {
		const y = pad.t + ph - t / top * ph;
		ctx.strokeStyle = t === 0 ? T.axis : T.grid;
		line(ctx, pad.l, y + 0.5, pad.l + pw, y + 0.5);
		ctx.fillStyle = T.muted;
		ctx.fillText(fmt(t, 3), pad.l - 6, y);
	}
	const bw = pw / picked.length;
	const gap = bw >= 6 ? 2 : bw >= 3 ? 1 : 0;
	const hv = hover[cv.id];
	picked.forEach((i, k) => {
		const bh = Math.max(1, w[i] / top * ph);
		const x = pad.l + k * bw + gap / 2, y = pad.t + ph - bh;
		ctx.fillStyle = signs[i] > 0 ? T.a : T.b;
		ctx.globalAlpha = hv === null || hv === undefined || hv === k ? 1 : 0.55;
		ctx.beginPath();
		if (ctx.roundRect && bw - gap >= 6) ctx.roundRect(x, y, bw - gap, bh, [3, 3, 0, 0]);
		else ctx.rect(x, y, Math.max(0.6, bw - gap), bh);
		ctx.fill();
	});
	ctx.globalAlpha = 1;
	cv._geom = {
		render: renderPartition,
		pick: (px) => {
			const k = Math.floor((px - pad.l) / bw);
			return k >= 0 && k < picked.length ? k : null;
		},
		tip: (k) => {
			const i = picked[k];
			return `<b>Number ${i}</b><div>a = ${fmt(w[i], 5)}</div><div><span class="sw" style="background:${signs[i] > 0 ? T.a : T.b}"></span>set ${signs[i] > 0 ? 'A (s = +1)' : 'B (s = −1)'}</div>`;
		},
	};

	// totals of the two sets
	let sa = 0, sb = 0;
	w.forEach((v, i) => { if (signs[i] > 0) sa += v; else sb += v; });
	let resid = Math.abs(sa - sb);
	if (setsMode === 'final') { sa = part.sum_a; sb = part.sum_b; resid = part.residual; }
	else if (setsMode === 'local') resid = part.residual_local;
	else if (data.residuals) resid = data.residuals.rounded[layer];
	$('problemNote').innerHTML = totalsHtml(sa, sb, resid, part);
}

function totalsHtml(sa, sb, resid, part) {
	const big = Math.max(sa, sb, 1e-12);
	const opt = part.optimum;
	return `
		<div class="totals">
			<span>Set A</span><span class="bar" style="width:${sa / big * 100}%;background:var(--set-a)"></span><span class="num">${fmt(sa, 6)}</span>
			<span>Set B</span><span class="bar" style="width:${sb / big * 100}%;background:var(--set-b)"></span><span class="num">${fmt(sb, 6)}</span>
		</div>
		<p class="verdict">Imbalance |A − B| = <b>${fmt(resid, 4)}</b> for this split. Karmarkar–Karp finds ${fmt(part.karmarkar_karp, 4)}${
			opt === null || opt === undefined ? '' : `, and the best possible is ${opt === 0 ? '0 (a perfect split)' : fmt(opt, 4)}`}.</p>`;
}

function renderSK() {
	const cv = $('problemCanvas'), j = data.problem.j, n = data.problem.n;
	problemCard('The coupling matrix J');
	if (!j) {
		cv.hidden = true;
		$('problemCaption').textContent = '';
		$('problemNote').innerHTML = `<p class="caption">J has ${n}×${n} = ${(n * n).toLocaleString()} entries; the picture is drawn for N ≤ 200. Run with a smaller N to see it.</p>`;
		return;
	}
	const signs = currentSigns() || finalAssignment(false);
	const order = Array.from({ length: n }, (_, i) => i).sort((a, b) => (signs[b] - signs[a]) || (a - b));
	const nUp = signs.filter((s) => s > 0).length;
	const q = j.q;
	let qmax = 1;
	for (const v of q) if (Math.abs(v) > qmax) qmax = Math.abs(v);

	// one pixel per entry, scaled up without smoothing
	const img = new ImageData(n, n);
	const mid = rgb(T.mid);
	let sat = 0, tot = 0;
	for (let r = 0; r < n; r++) {
		for (let c = 0; c < n; c++) {
			const i = order[r], k = order[c], v = q[i * n + k];
			let col = mid;
			if (skView === 'j') {
				if (v) col = mixRgb(T.mid, v > 0 ? T.pos : T.neg, Math.min(1, 0.35 + 0.65 * Math.abs(v) / qmax));
			} else if (v && i !== k) {
				const e = v * signs[i] * signs[k];
				col = mixRgb(T.mid, e < 0 ? T.neg : T.pos, Math.min(1, 0.35 + 0.65 * Math.abs(v) / qmax));
				if (r < c) { tot++; if (e < 0) sat++; }
			}
			const o = 4 * (r * n + c);
			img.data[o] = col[0]; img.data[o + 1] = col[1]; img.data[o + 2] = col[2]; img.data[o + 3] = 255;
		}
	}
	if (skView === 'j') {
		for (let r = 0; r < n; r++) for (let c = r + 1; c < n; c++) {
			const v = q[order[r] * n + order[c]];
			if (v) { tot++; if (v * signs[order[r]] * signs[order[c]] < 0) sat++; }
		}
	}
	const off = document.createElement('canvas');
	off.width = n; off.height = n;
	off.getContext('2d').putImageData(img, 0, 0);

	const side = Math.min(cv.clientWidth, 440);
	cv.style.height = side + 'px';
	const { ctx, w: W } = setup(cv);
	const size = side - 2, x0 = Math.floor((W - size) / 2), y0 = 1, cell = size / n;
	ctx.imageSmoothingEnabled = false;
	ctx.drawImage(off, x0, y0, size, size);
	// boundary between the two groups
	ctx.strokeStyle = T.text;
	ctx.lineWidth = 1.5;
	if (nUp > 0 && nUp < n) {
		line(ctx, x0 + nUp * cell, y0, x0 + nUp * cell, y0 + size);
		line(ctx, x0, y0 + nUp * cell, x0 + size, y0 + nUp * cell);
	}
	const hv = hover[cv.id];
	if (hv !== null && hv !== undefined) {
		ctx.strokeStyle = T.text;
		ctx.lineWidth = 1;
		ctx.strokeRect(x0 + (hv % n) * cell, y0 + Math.floor(hv / n) * cell, Math.max(cell, 2), Math.max(cell, 2));
	}

	const which = setsMode === 'layer' ? `the signs at layer ${layer}` : setsMode === 'local' ? 'the signs after the local search' : 'the final signs';
	$('problemCaption').innerHTML = skView === 'j'
		? `Rows and columns are spins, sorted by ${which}: s = +1 first (${nUp}), then s = −1 (${n - nUp}). Red: J<sub>ij</sub> &gt; 0, the pair prefers opposite signs; blue: J<sub>ij</sub> &lt; 0, it prefers equal signs.`
		: `Each pair's energy J<sub>ij</sub> s<sub>i</sub> s<sub>j</sub> for ${which}. Blue: the coupling is satisfied (lowers E); red: frustrated (raises E).`;
	$('problemNote').innerHTML = `<p class="verdict"><b>${(100 * sat / Math.max(1, tot)).toFixed(1)}%</b> of the ${tot.toLocaleString()} couplings are satisfied by this split (a random split satisfies about 50%).</p>`;
	cv._geom = {
		render: renderSK,
		pick: (px, py) => {
			const c = Math.floor((px - x0) / cell), r = Math.floor((py - y0) / cell);
			return c >= 0 && c < n && r >= 0 && r < n ? r * n + c : null;
		},
		tip: (idx) => {
			const i = order[Math.floor(idx / n)], k = order[idx % n], v = q[i * n + k];
			return `<b>Spins ${i} and ${k}</b><div>J = ${fmt(v * j.scale, 4)}</div><div>sᵢ sⱼ = ${signs[i] * signs[k] > 0 ? '+1' : '−1'}</div>` +
				(i === k ? '' : `<div>${v * signs[i] * signs[k] < 0 ? 'satisfied' : v ? 'frustrated' : 'no coupling'}</div>`);
		},
	};
}

function splitWords() {
	return setsMode === 'layer' ? `the signs at layer ${layer}` : setsMode === 'local' ? 'the final signs after the local search' : 'the final signs';
}

function renderMaxCut() {
	const cv = $('problemCanvas'), E = data.problem.edges, n = data.problem.n, mc = data.result.maxcut;
	problemCard('The graph and the cut');
	const signs = currentSigns() || finalAssignment(false);
	if (!E || !signs) {
		cv.hidden = true;
		$('problemCaption').textContent = '';
		$('problemNote').innerHTML = `<p class="caption">The graph is drawn for N ≤ 500. This run cut <b>${pct(mc.fraction)}</b> of its ${mc.edges.toLocaleString()} edges (mean field), ${pct(mc.local_fraction)} after the local search.</p>`;
		return;
	}
	if (!data._adj) {
		data._adj = Array.from({ length: n }, () => []);
		for (const [a, b] of E) { data._adj[a].push(b); data._adj[b].push(a); }
	}
	const W0 = cv.clientWidth;
	const H0 = Math.round(Math.max(260, Math.min(430, W0 * 0.62)));
	cv.style.height = H0 + 'px';
	const { ctx, w, h } = setup(cv);
	const sides = [[], []];
	for (let i = 0; i < n; i++) sides[signs[i] > 0 ? 0 : 1].push(i);
	const R = Math.max(30, Math.min(w / 4 - 18, h / 2 - 30));
	const cy = h / 2 + 10, cxs = [w / 4, 3 * w / 4];
	const pos = new Array(n);
	sides.forEach((list, s) => list.forEach((i, k) => {
		const a = 2 * Math.PI * k / Math.max(1, list.length) - Math.PI / 2;
		pos[i] = [cxs[s] + R * Math.cos(a), cy + R * Math.sin(a)];
	}));
	let cut = 0;
	for (const [a, b] of E) if (signs[a] !== signs[b]) cut++;
	const hv = hover[cv.id];
	const near = hv !== null && hv !== undefined ? new Set([hv, ...data._adj[hv]]) : null;

	// cut edges faint underneath, uncut edges (the ones that cost energy) on top
	for (const wantCut of [true, false]) {
		ctx.strokeStyle = wantCut ? T.muted : T.pos;
		ctx.lineWidth = wantCut ? 1 : 1.6;
		for (const [a, b] of E) {
			if ((signs[a] !== signs[b]) !== wantCut) continue;
			const touch = near && (a === hv || b === hv);
			ctx.globalAlpha = near ? (touch ? 1 : 0.12) : wantCut ? 0.35 : 0.9;
			line(ctx, pos[a][0], pos[a][1], pos[b][0], pos[b][1]);
		}
	}
	ctx.globalAlpha = 1;
	const r = n <= 60 ? 5 : n <= 200 ? 3.6 : 2.6;
	for (let i = 0; i < n; i++) {
		ctx.globalAlpha = near && !near.has(i) ? 0.3 : 1;
		dot(ctx, pos[i][0], pos[i][1], i === hv ? r + 2 : r, signs[i] > 0 ? T.a : T.b);
	}
	ctx.globalAlpha = 1;
	ctx.textAlign = 'center';
	ctx.textBaseline = 'top';
	ctx.font = '12px ' + getComputedStyle(document.body).fontFamily;
	haloText(ctx, `Side A (s = +1): ${sides[0].length}`, cxs[0], 4, T.a);
	haloText(ctx, `Side B (s = −1): ${sides[1].length}`, cxs[1], 4, T.b);

	$('problemCaption').innerHTML = `Vertices sorted onto two circles by ${splitWords()}. Edges between the circles are cut; red edges join two vertices on the same side, and each one costs energy.`;
	$('problemLegend').innerHTML = key(T.muted, 'cut edge', 'line') + key(T.pos, 'uncut edge (same side)', 'line') + '<span>hover a vertex for its edges</span>';
	$('problemNote').innerHTML = `<p class="verdict"><b>${cut}</b> of ${E.length} edges cut (<b>${pct(cut / E.length)}</b>) by this split. QAOA cuts 69.24% with one layer and 75.59% with two; mean field plus local search reaches ${pct(mc.local_fraction)}.</p>`;
	cv._geom = {
		render: renderMaxCut,
		pick: (px, py) => {
			let best = null, bd = 10;
			for (let i = 0; i < n; i++) {
				const d = Math.hypot(pos[i][0] - px, pos[i][1] - py);
				if (d < bd) { bd = d; best = i; }
			}
			return best;
		},
		tip: (i) => {
			const nb = data._adj[i], c = nb.filter((j) => signs[j] !== signs[i]).length;
			return `<b>Vertex ${i}</b><div><span class="sw" style="background:${signs[i] > 0 ? T.a : T.b}"></span>side ${signs[i] > 0 ? 'A (s = +1)' : 'B (s = −1)'}</div><div>${c} of its ${nb.length} edges cut</div>`;
		},
	};
}

function renderChain() {
	const cv = $('problemCanvas'), n = data.problem.n, ch = data.result.chain;
	problemCard('The chain');
	let z = null;
	if (setsMode === 'layer' && allShown()) z = layerZ(layer);
	else { const a = finalAssignment(setsMode === 'local'); if (a) z = a; }
	if (!z) {
		cv.hidden = true;
		$('problemCaption').textContent = '';
		$('problemNote').innerHTML = `<p class="caption">The sites are drawn for N ≤ 20,000. This run left <b>${ch.defects.toLocaleString()}</b> defects (mean field)${ch.defects_local === null || ch.defects_local === undefined ? '' : `, ${ch.defects_local.toLocaleString()} after the local search`}.</p>`;
	} else {
		const MAX = 2000;
		const m = Math.min(z.length, MAX);
		const W0 = cv.clientWidth;
		const rows = Math.max(1, Math.ceil(m * 11 / W0));
		const perRow = Math.ceil(m / rows);
		const cell = Math.min(24, (W0 - 2) / perRow);
		const rowH = cell + 14;
		cv.style.height = Math.round(rows * rowH + 4) + 'px';
		const { ctx } = setup(cv);
		const x0 = 1, sgn = (v) => (v < 0 ? -1 : 1);
		let defects = 0;
		for (let i = 0; i + 1 < z.length; i++) if (sgn(z[i]) === sgn(z[i + 1])) defects++;
		const hv = hover[cv.id];
		for (let i = 0; i < m; i++) {
			const r = Math.floor(i / perRow), c = i % perRow;
			const x = x0 + c * cell, y = 2 + r * rowH;
			const v = Math.max(-1, Math.min(1, z[i]));
			ctx.fillStyle = css(mixRgb(T.mid, v >= 0 ? T.a : T.b, Math.min(1, 0.15 + 0.85 * Math.abs(v))));
			ctx.fillRect(x + 0.5, y, Math.max(1, cell - 1), cell);
			// which of the two alternating patterns the site follows
			const d = sgn(z[i]) * (i % 2 ? -1 : 1);
			ctx.fillStyle = d > 0 ? T.c3 : T.c7;
			ctx.globalAlpha = 0.25 + 0.75 * Math.abs(v);
			ctx.fillRect(x + 0.5, y + cell + 2, Math.max(1, cell - 1), 4);
			ctx.globalAlpha = 1;
			if (i + 1 < z.length && sgn(z[i]) === sgn(z[i + 1])) {
				const bx = c === perRow - 1 ? x + cell - 1 : x + cell;
				ctx.fillStyle = T.pos;
				ctx.beginPath();
				ctx.moveTo(bx - 4, y + cell + 12);
				ctx.lineTo(bx + 4, y + cell + 12);
				ctx.lineTo(bx, y + cell + 6);
				ctx.fill();
			}
			if (i === hv) {
				ctx.strokeStyle = T.text;
				ctx.lineWidth = 1.5;
				ctx.strokeRect(x + 0.5, y, cell - 1, cell);
			}
		}
		const what = setsMode === 'layer' ? `z at layer ${layer}` : setsMode === 'local' ? 'the final signs after the local search' : 'the final signs';
		$('problemCaption').innerHTML = `Each cell is one site, colored by ${what}${m < z.length ? ` (the first ${m.toLocaleString()} of ${z.length.toLocaleString()} sites)` : ''}. The thin band shows which of the two alternating patterns, ↑↓↑↓ or ↓↑↓↑, each site follows: where it changes color there is a domain wall. Red marks: two neighbours alike (a defect).`;
		$('problemLegend').innerHTML = key(T.a, 'z &gt; 0', '') + key(T.b, 'z &lt; 0', '') + key(T.c3, 'follows ↑↓↑↓', '') + key(T.c7, 'follows ↓↑↓↑', '') + key(T.pos, 'defect', '');
		$('problemNote').innerHTML = `<p class="verdict"><b>${defects.toLocaleString()}</b> defect${defects === 1 ? '' : 's'} in this state; the ground state has none. At the end: ${ch.defects} with mean field${ch.defects_local === null || ch.defects_local === undefined ? '' : `, ${ch.defects_local} after the local search`} (single flips can't move a domain wall, so local search rarely removes them). Energy ${fmt(data.result.rounded_energy, 6)} against ${fmt(ch.ground_energy, 6)} for perfect alternation.</p>`;
		cv._geom = {
			render: renderChain,
			pick: (px, py) => {
				const r = Math.floor((py - 2) / rowH), c = Math.floor((px - x0) / cell);
				const i = r * perRow + c;
				return r >= 0 && c >= 0 && c < perRow && i < m ? i : null;
			},
			tip: (i) => `<b>Site ${i}</b><div>${setsMode === 'layer' ? 'z = ' + fmt(z[i], 4) : 's = ' + (z[i] > 0 ? '+1' : '−1')}</div><div>follows ${sgn(z[i]) * (i % 2 ? -1 : 1) > 0 ? '↑↓↑↓' : '↓↑↓↑'}</div>`,
		};
	}

	// the couplings against distance
	const f = data.problem.couplings;
	if (f && f.length > 2) {
		const cv2 = $('problemCanvas2');
		cv2.hidden = false;
		$('problemCaption2').hidden = false;
		$('problemCaption2').innerHTML = `Coupling J = 1/r<sup>α</sup> against distance r (α = ${fmt(data.problem.alpha, 3)}), log scale.`;
		lineChart(cv2, {
			series: [{ name: 'J(r)', values: f.slice(1), color: T.c7, x0: 1 }],
			xmin: 1, logY: true, floor: 1e-300, xName: 'Distance', render: renderChain,
		});
	}
}

function renderHopfield() {
	const cv = $('problemCanvas'), pr = data.problem, n = pr.n, hp = data.result.hopfield;
	problemCard('Stored patterns and the state');
	const P = pr.patterns, cue = pr.cue;
	let state = null, stateTitle = '';
	if (setsMode === 'layer' && allShown()) { state = layerZ(layer); stateTitle = `Layer ${layer}: z`; }
	else {
		state = finalAssignment(setsMode === 'local');
		stateTitle = setsMode === 'local' ? 'Final + local search' : 'Final: sign of z';
	}
	const overlap = (x, s) => { let t = 0; for (let i = 0; i < n; i++) t += x[i] * s[i]; return t / n; };
	const best = hp.best_pattern, inv = hp.best_overlap < 0;
	const panels = [];
	if (P) panels.push({ title: 'Stored pattern 1', vals: P[0] });
	if (cue && P) panels.push({ title: `Cue (${Math.round(100 * (pr.noise ?? 0))}% flipped)`, vals: cue, sub: `overlap ${fmt(overlap(P[0], cue), 3)}` });
	else if (cue) panels.push({ title: 'Cue', vals: cue });
	if (state) panels.push({ title: stateTitle, vals: state, sub: P ? `overlap with 1: ${fmt(overlap(P[0], state), 3)}` : '' });
	if (P && P[best - 1]) {
		const b = inv ? P[best - 1].map((v) => -v) : P[best - 1];
		panels.push({ title: `Best match: pattern ${best}${inv ? ', inverted' : ''}`, vals: b, sub: `final overlap ${fmt(Math.abs(hp.best_overlap), 3)}` });
	}
	if (!panels.length || n > 40000) {
		cv.hidden = true;
		$('problemCaption').textContent = 'The patterns are drawn when N × patterns ≤ 100,000.';
	} else {
		const rows0 = pr.shape ? pr.shape[0] : Math.ceil(n / Math.ceil(Math.sqrt(n)));
		const cols0 = pr.shape ? pr.shape[1] : Math.ceil(Math.sqrt(n));
		const W0 = cv.clientWidth;
		const cols = W0 >= 520 ? panels.length : Math.min(2, panels.length);
		const gap = 14;
		const size = Math.min(170, (W0 - (cols - 1) * gap) / cols);
		const blockH = size + 44;
		const prow = Math.ceil(panels.length / cols);
		cv.style.height = Math.round(prow * blockH) + 'px';
		const { ctx, w } = setup(cv);
		const totalW = cols * size + (cols - 1) * gap, left = (w - totalW) / 2;
		panels.forEach((p, k) => {
			const r = Math.floor(k / cols), c = k % cols;
			const x = left + c * (size + gap), y = r * blockH + 18;
			ctx.textAlign = 'center';
			ctx.textBaseline = 'bottom';
			ctx.fillStyle = T.text2;
			ctx.font = '12px ' + getComputedStyle(document.body).fontFamily;
			ctx.fillText(p.title, x + size / 2, y - 4, size + gap);
			const cell = size / Math.max(rows0, cols0);
			drawPixels(ctx, p.vals, rows0, cols0, x + (size - cols0 * cell) / 2, y, size);
			if (p.sub) {
				ctx.textBaseline = 'top';
				ctx.fillStyle = T.muted;
				ctx.font = '11px ' + getComputedStyle(document.body).fontFamily;
				ctx.fillText(p.sub, x + size / 2, y + rows0 * cell + 5, size + gap);
			}
		});
		$('problemCaption').innerHTML = `Each image is ${n} spins${pr.shape ? '' : ' laid out in rows'}: dark +1, light −1${setsMode === 'layer' ? ' (in between: z near 0)' : ''}. The patterns are random, so a match shows as the same speckle.`;
	}

	// overlaps with every stored pattern
	const ov = hp.overlaps;
	if (ov && ov.length) {
		const cv2 = $('problemCanvas2');
		cv2.hidden = false;
		$('problemCaption2').hidden = false;
		$('problemCaption2').innerHTML = `Overlap of the final mean-field state with each stored pattern: ±1 means equal or inverted, near 0 unrelated.${ov.length < hp.patterns ? ` First ${ov.length} of ${hp.patterns}.` : ''}`;
		renderOverlapBars(cv2, ov, best);
	}
	$('problemNote').innerHTML = loadGauge(hp) + `<p class="verdict">Final state: overlap <b>${fmt(Math.abs(hp.best_overlap), 3)}</b> with pattern ${best}${inv ? ' (inverted)' : ''}${Math.abs(hp.best_overlap) > 0.9 ? ', a stored memory' : ', a mixture rather than a stored memory'}. Classic recall from the cue reaches overlap ${fmt(hp.overlap_recall, 3)} with pattern 1.</p>`;
}

function renderOverlapBars(cv, ov, best) {
	const { ctx, w, h } = setup(cv);
	const pad = { l: 40, r: 6, t: 8, b: 20 };
	const pw = w - pad.l - pad.r, ph = h - pad.t - pad.b;
	const Y = (v) => pad.t + (1 - v) / 2 * ph;
	ctx.textAlign = 'right';
	ctx.textBaseline = 'middle';
	ctx.lineWidth = 1;
	for (const t of [-1, -0.5, 0, 0.5, 1]) {
		ctx.strokeStyle = t === 0 ? T.axis : T.grid;
		line(ctx, pad.l, Y(t) + 0.5, pad.l + pw, Y(t) + 0.5);
		ctx.fillStyle = T.muted;
		ctx.fillText(String(t), pad.l - 6, Y(t));
	}
	const bw = pw / ov.length, gap = bw >= 6 ? 2 : 0, hv = hover[cv.id];
	ov.forEach((v, k) => {
		ctx.fillStyle = k === 0 ? T.a : k === best - 1 ? T.b : T.muted;
		ctx.globalAlpha = hv === null || hv === undefined || hv === k ? 1 : 0.55;
		const y0 = Y(0), y1 = Y(v);
		ctx.fillRect(pad.l + k * bw + gap / 2, Math.min(y0, y1), Math.max(0.8, bw - gap), Math.max(1, Math.abs(y1 - y0)));
	});
	ctx.globalAlpha = 1;
	ctx.textAlign = 'center';
	ctx.textBaseline = 'top';
	ctx.fillStyle = T.muted;
	for (const t of niceTicks(1, ov.length, Math.min(8, ov.length))) {
		if (t !== Math.round(t)) continue;
		ctx.fillText(String(t), pad.l + (t - 0.5) * bw, pad.t + ph + 4);
	}
	cv._geom = {
		render: () => renderOverlapBars(cv, ov, best),
		pick: (px) => { const k = Math.floor((px - pad.l) / bw); return k >= 0 && k < ov.length ? k : null; },
		tip: (k) => `<b>Pattern ${k + 1}</b><div>overlap ${fmt(ov[k], 3)}</div>${k === 0 ? '<div>the cue\'s pattern</div>' : ''}${k === best - 1 ? '<div>best match</div>' : ''}`,
	};
}

function loadGauge(hp) {
	const top = Math.max(0.2, hp.load * 1.15);
	const at = (v) => (100 * Math.min(v, top) / top).toFixed(2) + '%';
	const g = hp.global_limit ?? 0.051, c = hp.capacity ?? 0.138;
	const verdict = hp.load < g ? 'below 0.05: the stored patterns are the lowest-energy states, so annealing can find them'
		: hp.load < c ? 'between 0.05 and 0.138: recall from a cue works, but annealing ends in mixtures'
			: 'above 0.138: even recall from a cue fails';
	return `<div class="gauge">
		<div class="gtrack">
			<span class="gzone ok" style="left:0;width:${at(g)}"></span>
			<span class="gzone mid" style="left:${at(g)};width:calc(${at(c)} - ${at(g)})"></span>
			<span class="gzone bad" style="left:${at(c)};right:0"></span>
			<span class="gmark" style="left:${at(hp.load)}"></span>
		</div>
		<div class="gticks"><span style="left:0">0</span><span style="left:${at(g)}">0.051</span><span style="left:${at(c)}">0.138</span></div>
		<p class="verdict">Load P/N = <b>${fmt(hp.load, 3)}</b>: ${verdict}.</p>
	</div>`;
}


// ---------------------------------------------------------------
// Bloch sphere
// ---------------------------------------------------------------

let drag = null;

function project(p, cx, cy, R) {
	const cyaw = Math.cos(view.yaw), syaw = Math.sin(view.yaw);
	const x = p[0] * cyaw - p[1] * syaw, y = p[0] * syaw + p[1] * cyaw, z = p[2];
	const cp = Math.cos(view.pitch), sp = Math.sin(view.pitch);
	return { x: cx + R * y, y: cy - R * (-x * sp + z * cp), d: x * cp + z * sp };
}

function renderSphere() {
	const cv = $('sphere');
	const { ctx, w, h } = setup(cv);
	const R = Math.min(w, h) / 2 - 24, cx = w / 2, cy = h / 2;
	const P = (p) => project(p, cx, cy, R);

	ctx.strokeStyle = T.axis;
	ctx.lineWidth = 1;
	ctx.beginPath();
	ctx.arc(cx, cy, R, 0, 2 * Math.PI);
	ctx.stroke();

	// equator and two meridians, back halves faint
	const curve = (f) => {
		let prev = P(f(0));
		for (let k = 1; k <= 96; k++) {
			const cur = P(f(k / 96 * 2 * Math.PI));
			ctx.globalAlpha = (prev.d + cur.d) / 2 < 0 ? 0.35 : 1;
			line(ctx, prev.x, prev.y, cur.x, cur.y);
			prev = cur;
		}
		ctx.globalAlpha = 1;
	};
	ctx.strokeStyle = T.grid;
	curve((t) => [Math.cos(t), 0, Math.sin(t)]);
	curve((t) => [0, Math.cos(t), Math.sin(t)]);
	ctx.strokeStyle = T.axis;
	curve((t) => [Math.cos(t), Math.sin(t), 0]);

	// axes
	ctx.textAlign = 'center';
	ctx.textBaseline = 'middle';
	ctx.font = '12px ' + getComputedStyle(document.body).fontFamily;
	for (const [e, name] of [[[1, 0, 0], 'x'], [[0, 1, 0], 'y'], [[0, 0, 1], 'z']]) {
		const a = P(e.map((v) => -v)), b = P(e), c = P(e.map((v) => v * 1.14));
		ctx.strokeStyle = T.muted;
		ctx.setLineDash([3, 3]);
		line(ctx, a.x, a.y, cx, cy);
		ctx.setLineDash([]);
		line(ctx, cx, cy, b.x, b.y);
		haloText(ctx, name, c.x, c.y, T.text2);
	}

	if (!data) return;
	const sp = data.layers[layer].spins, m = data.shown.length;
	const fa = finalAssignment(false);
	const pts = [];
	for (let k = 0; k < m; k++) {
		const p = P([sp[3 * k], sp[3 * k + 1], sp[3 * k + 2]]);
		const s = fa ? fa[data.shown[k]] : (sp[3 * k + 2] < 0 ? -1 : 1);
		pts.push({ ...p, color: s > 0 ? T.a : T.b });
	}
	pts.sort((a, b) => a.d - b.d);
	const r = m <= 60 ? 5 : m <= 300 ? 3.5 : 2.5;
	ctx.lineWidth = 1;
	for (const p of pts) {
		ctx.globalAlpha = p.d < 0 ? 0.45 : 1;
		ctx.beginPath();
		ctx.arc(p.x, p.y, r, 0, 2 * Math.PI);
		ctx.fillStyle = p.color;
		ctx.fill();
		ctx.strokeStyle = T.surface;
		ctx.stroke();
	}
	ctx.globalAlpha = 1;
}

function sphereLegend() {
	const fa = finalAssignment(false);
	const kind = data.problem.kind;
	const word = kind === 'partition' ? ['set A', 'set B'] : kind === 'maxcut' ? ['side A', 'side B'] : ['s = +1', 's = −1'];
	$('sphereLegend').innerHTML = (fa
		? key(T.a, `ends with z &gt; 0 (${word[0]})`, 'dot') + key(T.b, `ends with z &lt; 0 (${word[1]})`, 'dot')
		: key(T.a, 'z &gt; 0 now', 'dot') + key(T.b, 'z &lt; 0 now', 'dot')) +
		(data.shown.length < data.problem.n ? `<span>${data.shown.length} of ${data.problem.n} spins drawn</span>` : '') +
		'<span>faint: far side of the sphere</span>';
}

function attachSphere() {
	const cv = $('sphere');
	cv.addEventListener('pointerdown', (e) => {
		drag = { x: e.clientX, y: e.clientY };
		cv.setPointerCapture(e.pointerId);
	});
	cv.addEventListener('pointermove', (e) => {
		if (!drag) return;
		view.yaw += (e.clientX - drag.x) * 0.009;
		view.pitch = Math.max(-1.5, Math.min(1.5, view.pitch + (e.clientY - drag.y) * 0.009));
		drag = { x: e.clientX, y: e.clientY };
		renderSphere();
	});
	const end = () => { drag = null; };
	cv.addEventListener('pointerup', end);
	cv.addEventListener('pointercancel', end);
}


// ---------------------------------------------------------------
// charts over layers
// ---------------------------------------------------------------

function renderEnergy() {
	const L = data.layers, hw = data.params.engine === 'hw';
	const series = [{ name: hw ? 'mean field, bit-exact hw' : 'mean field', values: L.map((l) => l.energy), color: T.a }];
	if (data.compare) series.push({ name: 'mean field, floating point', values: data.compare.energy, color: T.c3, under: true });
	series.push({ name: 'rounded to ±1', values: L.map((l) => l.rounded_energy), color: T.b });
	lineChart($('energy'), {
		series,
		hlines: data.baselines.map((b) => ({ label: b.name, value: b.energy })),
		marker: layer,
		render: renderEnergy,
		click: setLayer,
	});
	const term = {
		partition: 'a<sub>i</sub>a<sub>j</sub>z<sub>i</sub>z<sub>j</sub>',
		chain: 'z<sub>i</sub>z<sub>j</sub>/|i − j|<sup>α</sup>',
	}[data.problem.kind] || 'J<sub>ij</sub>z<sub>i</sub>z<sub>j</sub>';
	$('energyCaption').innerHTML = `Mean-field energy Σ<sub>i&lt;j</sub> ${term} of the spin vectors, and the energy after rounding each z to ±1.`;
	$('energyLegend').innerHTML = series.map((s) => key(s.color, s.name, 'line')).join('') +
		key(null, 'classical methods / limits', 'dash');
}

// the problem's own measure per layer; older files only have the partition residuals
function metricsOf(d) {
	if (d.metrics) return d.metrics;
	if (!d.residuals) return null;
	const p = d.result.partition;
	return {
		title: 'Imbalance |sum A − sum B|', log: true,
		mean_field: { name: 'mean field |Σ aᵢ zᵢ|', values: d.residuals.mean_field },
		rounded: { name: 'rounded |Σ aᵢ sᵢ|', values: d.residuals.rounded },
		refs: [['mean field + local search', p.residual_local], ['Karmarkar–Karp', p.karmarkar_karp], ['greedy', p.greedy], ['optimum (exact)', p.optimum]]
			.map(([name, value]) => ({ name, value })),
	};
}

function renderMetric() {
	const m = metricsOf(data), card = $('metricCard');
	if (!m) { card.hidden = true; return; }
	card.hidden = false;
	const kind = data.problem.kind;
	let floor = 1e-12;
	if (m.log) {
		let lo = Infinity;
		for (const v of [...m.mean_field.values, ...m.rounded.values, ...m.refs.map((r) => r.value)]) if (v > 0 && v < lo) lo = v;
		floor = kind === 'partition' ? 1 / 32768 : isFinite(lo) ? lo / 2 : 1e-12;
	}
	const clip = (v) => (m.log ? Math.max(v, floor) : v);
	const series = [
		{ name: m.mean_field.name, values: m.mean_field.values.map(clip), color: T.a },
		{ name: m.rounded.name, values: m.rounded.values.map(clip), color: T.b },
	];
	lineChart($('metric'), {
		series,
		hlines: m.refs.map((r) => ({ label: r.name, value: r.value })),
		logY: m.log, floor, marker: layer, render: renderMetric, click: setLayer,
	});
	$('metricTitle').textContent = m.title + ' per layer';
	const better = { partition: 'Lower is better', sk: 'Lower is better', 'sk-gauss': 'Lower is better' }[kind] || 'Higher is better';
	$('metricCaption').textContent = `${m.log ? 'Log scale. ' : ''}${better}. The mean-field line uses the spin vectors, the other rounds each z to ±1; dashed: classical methods and limits.`;
	const zeroRef = m.refs.find((r) => r.value === 0 && m.log);
	$('metricLegend').innerHTML = series.map((s) => key(s.color, s.name, 'line')).join('') + key(null, 'classical methods / limits', 'dash') +
		(zeroRef ? `<span>${zeroRef.name}: 0, which a log axis can't show</span>` : '');
}

function renderSchedule() {
	const s = data.schedule;
	const series = [
		{ name: 'γ (cost)', values: s.gamma, color: T.c3, x0: 1 },
		{ name: 'β (mixer)', values: s.beta, color: T.c7, x0: 1 },
	];
	lineChart($('schedule'), {
		series, xmin: 1, marker: layer >= 1 ? layer : null, render: renderSchedule, click: setLayer,
	});
	$('scheduleCaption').textContent = data.params.schedule === 'ramp'
		? 'Cost step strength γ rises and mixer strength β falls, like a slow anneal.'
		: 'The hardware\'s angle_lut values: a test pattern, not an optimizing schedule.';
	$('scheduleLegend').innerHTML = series.map((x) => key(x.color, x.name, 'line')).join('') +
		'<span>layer k turns spins by 2γ<sub>k</sub>h<sub>i</sub>/ĥ about z, then by 2β<sub>k</sub> about x</span>';
}

function renderHist() {
	const cv = $('hist');
	const { ctx, w, h } = setup(cv);
	const counts = data.layers[layer].z_hist, bins = counts.length;
	if (!data._histMax) data._histMax = Math.max(1, ...data.layers.map((l) => Math.max(...l.z_hist)));
	const top = Math.max(...niceTicks(0, data._histMax, 4));
	const pad = { l: 44, r: 8, t: 8, b: 24 };
	const pw = w - pad.l - pad.r, ph = h - pad.t - pad.b;
	ctx.textAlign = 'right';
	ctx.textBaseline = 'middle';
	ctx.lineWidth = 1;
	for (const t of niceTicks(0, top, 4)) {
		const y = pad.t + ph - t / top * ph;
		ctx.strokeStyle = t === 0 ? T.axis : T.grid;
		line(ctx, pad.l, y + 0.5, pad.l + pw, y + 0.5);
		ctx.fillStyle = T.muted;
		ctx.fillText(t.toLocaleString(), pad.l - 6, y);
	}
	ctx.textAlign = 'center';
	ctx.textBaseline = 'top';
	for (const t of [-1, -0.5, 0, 0.5, 1]) ctx.fillText(String(t), pad.l + (t + 1) / 2 * pw, pad.t + ph + 6);
	const bw = pw / bins, gap = bw >= 6 ? 2 : 1, hv = hover[cv.id];
	counts.forEach((c, k) => {
		if (!c) return;
		const bh = Math.max(1.5, c / top * ph), x = pad.l + k * bw + gap / 2;
		ctx.fillStyle = k < bins / 2 ? T.b : T.a;
		ctx.globalAlpha = hv === null || hv === undefined || hv === k ? 1 : 0.55;
		ctx.beginPath();
		if (ctx.roundRect && bw - gap >= 5) ctx.roundRect(x, pad.t + ph - bh, bw - gap, bh, [3, 3, 0, 0]);
		else ctx.rect(x, pad.t + ph - bh, bw - gap, bh);
		ctx.fill();
	});
	ctx.globalAlpha = 1;
	cv._geom = {
		render: renderHist,
		pick: (px) => {
			const k = Math.floor((px - pad.l) / bw);
			return k >= 0 && k < bins ? k : null;
		},
		tip: (k) => `<b>z from ${(-1 + 2 * k / bins).toFixed(2)} to ${(-1 + 2 * (k + 1) / bins).toFixed(2)}</b><div>${counts[k].toLocaleString()} spins</div>`,
	};
}


// ---------------------------------------------------------------
// results
// ---------------------------------------------------------------

function tile(k, v, s = '') {
	return `<div class="tile"><div class="k">${k}</div><div class="v">${v}</div><div class="s">${s}</div></div>`;
}

// a table of methods: the mean-field row first, then the baselines
function methodTable(metricName, metricFmt, extra) {
	const m = metricsOf(data);
	const last = m ? m.rounded.values[m.rounded.values.length - 1] : null;
	const rows = [{ name: 'mean field (rounded)', metric: last, energy: data.result.rounded_energy, mf: true }]
		.concat(data.baselines.map((b) => ({ name: b.name, metric: b.metric, energy: b.energy })));
	return `<table><thead><tr><th>Method</th><th class="num">${metricName}</th><th class="num">Energy</th>${extra ? `<th class="num">${extra.head}</th>` : ''}</tr></thead><tbody>` +
		rows.map((r) => `<tr class="${r.mf ? 'mf' : ''}"><td>${r.name}</td><td class="num">${metricFmt(r.metric)}</td><td class="num">${fmt(r.energy, 7)}</td>${extra ? `<td class="num">${extra.cell(r)}</td>` : ''}</tr>`).join('') +
		'</tbody></table>';
}

function renderResults() {
	const r = data.result, p = data.params, t = data.timing, n = data.problem.n;
	let html = '';
	if (r.partition) {
		const q = r.partition;
		html += '<div class="tiles">' +
			tile('Imbalance, mean field', fmt(q.residual, 4), 'signs of z after the last layer') +
			tile('+ local search', fmt(q.residual_local, 4), 'single moves that shrink it') +
			tile('Karmarkar–Karp', fmt(q.karmarkar_karp, 4), 'differencing heuristic') +
			tile('Optimum', q.optimum === null || q.optimum === undefined ? '—' : fmt(q.optimum, 4), q.optimum === null || q.optimum === undefined ? 'too large to compute exactly' : 'exact, by dynamic programming') +
			'</div>';
		const rows = [
			['mean field (rounded)', q.residual, r.rounded_energy, true],
			['mean field + local search', q.residual_local],
			['Karmarkar–Karp', q.karmarkar_karp],
			['greedy (largest first)', q.greedy],
		];
		if (q.optimum !== null && q.optimum !== undefined) rows.push(['optimum (exact)', q.optimum]);
		html += `<table><thead><tr><th>Method</th><th class="num">|sum A − sum B|</th><th class="num">Energy</th></tr></thead><tbody>` +
			rows.map(([name, res, e, mf]) => {
				const be = data.baselines.find((b) => b.name.startsWith(name.split(' (')[0]));
				const energy = e !== undefined ? e : be ? be.energy : null;
				return `<tr class="${mf ? 'mf' : ''}"><td>${name}</td><td class="num">${fmt(res, 6)}</td><td class="num">${fmt(energy, 8)}</td></tr>`;
			}).join('') +
			`</tbody></table><p class="note">Sum of all numbers ${fmt(q.total, 6)}; set A ${fmt(q.sum_a, 6)}, set B ${fmt(q.sum_b, 6)} (mean field). Final mean-field energy of the spin vectors: ${fmt(r.energy, 6)}.</p>`;
	} else if (r.sk) {
		const s = r.sk;
		html += '<div class="tiles">' +
			tile('Energy per spin, mean field', fmt(s.energy_per_spin, 4), 'signs of z after the last layer') +
			tile('+ local search', fmt(s.local_per_spin, 4), 'single spin flips downhill') +
			tile('Random start + local search', fmt(s.random_local_per_spin, 4), 'the same search from a random split') +
			tile('Best possible, N → ∞', fmt(s.parisi_per_spin, 4), 'Parisi; finite N sits a little above') +
			'</div>';
		const rows = [['mean field (rounded)', r.rounded_energy, true]].concat(data.baselines.map((b) => [b.name, b.energy]));
		html += `<table><thead><tr><th>Method</th><th class="num">Energy</th><th class="num">Per spin</th></tr></thead><tbody>` +
			rows.map(([name, e, mf]) => `<tr class="${mf ? 'mf' : ''}"><td>${name}</td><td class="num">${fmt(e, 6)}</td><td class="num">${fmt(e / n, 4)}</td></tr>`).join('') +
			`</tbody></table><p class="note">Final mean-field energy of the spin vectors: ${fmt(r.energy, 6)}.</p>`;
	} else if (r.maxcut) {
		const m = r.maxcut;
		html += '<div class="tiles">' +
			tile('Edges cut, mean field', pct(m.fraction), `${m.cut.toLocaleString()} of ${m.edges.toLocaleString()}`) +
			tile('+ local search', pct(m.local_fraction), 'single flips that cut more') +
			tile('Random start + local search', pct(m.random_local_fraction), 'the same search from a random split') +
			(m.degree === 3 ? tile('QAOA, 1 and 2 layers', '69.24% · 75.59%', 'on large random 3-regular graphs') : '') +
			'</div>';
		html += methodTable('Edges cut', (v) => pct(v), { head: 'Count', cell: (row) => (row.metric === null || row.metric === undefined ? '—' : Math.round(row.metric * m.edges).toLocaleString()) });
		html += `<p class="note">Energy = edges − 2 × cut. Final mean-field energy of the spin vectors: ${fmt(r.energy, 6)}.</p>`;
	} else if (r.chain) {
		const c = r.chain;
		html += '<div class="tiles">' +
			tile('Defects, mean field', c.defects.toLocaleString(), 'neighbours alike; 0 in the ground state') +
			tile('+ local search', c.defects_local === null || c.defects_local === undefined ? '—' : c.defects_local.toLocaleString(), 'single flips can\'t move walls') +
			tile('Neighbours alternating', pct(c.alternating, 1), `overall order ${fmt(c.order, 3)}`) +
			tile('Energy above ground', fmt(r.rounded_energy - c.ground_energy, 4), `ground state ${fmt(c.ground_energy, 6)}`) +
			'</div>';
		html += methodTable('Alternating', (v) => pct(v, 1), { head: 'Above ground', cell: (row) => fmt(row.energy - c.ground_energy, 4) });
		html += `<p class="note">Overall order |Σ (−1)<sup>i</sup> s<sub>i</sub>| / N is ${fmt(c.order, 3)}: domains of the two alternating patterns cancel, so a few defects lower it a lot. α = ${fmt(c.alpha, 3)}.</p>`;
	} else if (r.hopfield) {
		const hp = r.hopfield;
		html += '<div class="tiles">' +
			tile('Best match', `pattern ${hp.best_pattern}`, `overlap ${fmt(hp.best_overlap, 3)}${hp.best_overlap < 0 ? ' (inverted)' : ''}`) +
			tile('Overlap with pattern 1', fmt(hp.overlap, 3), `the cue starts at ${fmt(hp.overlap_cue, 3)}`) +
			tile('Classic recall from the cue', fmt(hp.overlap_recall, 3), 'local search from the cue') +
			tile('Load P/N', fmt(hp.load, 3), `limits ${hp.global_limit ?? 0.051} (anneal), ${hp.capacity} (recall)`) +
			'</div>';
		html += methodTable('Largest overlap', (v) => fmt(v, 3));
		html += `<p class="note">Overlap = Σ ξ<sub>i</sub> s<sub>i</sub> / N with a stored pattern ξ; the table gives the largest over all ${hp.patterns} patterns. Mean field + local search reaches overlap ${fmt(hp.overlap_local, 3)} with pattern 1.</p>`;
	}

	if (data.hardware) {
		const hw = data.hardware;
		const leds = hw.led.split('').map((b) => `<span class="led ${b === '1' ? 'on' : ''}"></span>`).join('');
		html += `<div class="subhead">What the FPGA design would show</div>
			<div class="hwrow">
				<span>Energy register: <b>${hw.energy_q15.toLocaleString()}</b> (Q1.15, = ${fmt(hw.energy_q15 / 32768, 6)})</span>
				<span>LEDs<span class="leds" title="${hw.led}">${leds}</span><code>${hw.led}</code></span>
				${data.compare ? `<span>Largest hw − float spin difference: <b>${fmt(data.compare.max_spin_diff, 3)}</b></span>` : ''}
			</div>`;
	}

	const f = t.fpga, ratio = t.cpu_seconds / f.seconds;
	html += `<div class="subhead">Time</div><div class="tiles">` +
		tile(`CPU: ${t.cpu}`, fmtTime(t.cpu_seconds), `${t.threads} thread${t.threads > 1 ? 's' : ''}, simulation only`) +
		tile('FPGA estimate', fmtTime(f.seconds), f.clocks ? `${f.clocks.toLocaleString()} clocks` : '') +
		tile(ratio >= 1 ? 'FPGA faster by' : 'CPU faster by', `${fmt(ratio >= 1 ? ratio : 1 / ratio, 3)}×`, 'for this run') +
		`</div><p class="note">${f.design}${f.note ? '. ' + f.note : ''}. Run: N = ${n}, ${p.layers} layers, ${p.schedule === 'ramp' ? 'dt ' + p.dt : 'angle_lut'}, engine ${p.engine}, seed ${p.seed}.</p>`;
	$('results').innerHTML = html;
}


// ---------------------------------------------------------------
// layers, playing, running
// ---------------------------------------------------------------

function layerLabel() {
	const L = data.layers.length - 1, e = data.layers[layer];
	$('layerlabel').textContent = `Layer ${layer} of ${L} · E = ${fmt(e.energy, 5)}`;
}

function setLayer(k) {
	if (!data) return;
	layer = Math.max(0, Math.min(data.layers.length - 1, k | 0));
	$('layer').value = layer;
	layerLabel();
	renderSphere();
	renderEnergy();
	renderMetric();
	renderSchedule();
	renderHist();
	if (setsMode === 'layer') renderProblem();
}

let playStart = 0, playFrom = 0;

function setPlaying(on) {
	playing = on;
	$('play').textContent = on ? '❚❚' : '▶';
	$('play').setAttribute('aria-label', on ? 'Pause' : 'Play');
	if (!on || !data) return;
	if (layer >= data.layers.length - 1) setLayer(0);
	playStart = performance.now();
	playFrom = layer;
	requestAnimationFrame(tick);
}

function tick(now) {
	if (!playing || !data) return;
	const L = data.layers.length - 1;
	const perLayer = Math.max(6, 6000 / Math.max(1, L));	// whole run in about 6 s
	const k = Math.min(L, playFrom + Math.floor((now - playStart) / perLayer));
	if (k !== layer) setLayer(k);
	if (k >= L) setPlaying(false); else requestAnimationFrame(tick);
}

function renderAll() {
	T = readTheme();
	if (mode === 'studies') { renderStudy(); return; }
	if (!data) { renderSphere(); return; }
	renderAbout();
	renderProblem();
	sphereLegend();
	renderSphere();
	renderEnergy();
	renderMetric();
	renderSchedule();
	renderHist();
	renderResults();
	layerLabel();
}

// fromFile: also set the controls to the run's parameters
function show(d, fromFile = false) {
	data = d;
	setPlaying(false);
	if (!allShown() && setsMode === 'layer') setsMode = 'final';
	if (setsMode === 'local' && !finalAssignment(true)) setsMode = 'final';
	$('layer').max = d.layers.length - 1;
	layer = d.layers.length - 1;		// start on the result; play runs from layer 0
	$('layer').value = layer;
	if (fromFile) fillControls(d.params);
	renderAll();
}

function setStatus(html, error = false) {
	const s = $('status');
	s.innerHTML = html;
	s.classList.toggle('error', error);
}

function query() {
	const problem = F('problem').value;
	const keys = ['problem', 'n', 'layers', 'dt', 'engine', 'schedule', 'init', 'tilt', 'seed']
		.concat(Object.keys(EXTRAS).filter((k) => EXTRAS[k] === problem));
	return keys.map((k) => `${k}=${encodeURIComponent(F(k).value)}`).join('&');
}

const fileNote = 'This page was opened as a file, so it can\'t run simulations. Start <code>viral-sim serve</code> and open <a href="http://127.0.0.1:8080/">http://127.0.0.1:8080/</a>, or load a JSON file.';

async function run() {
	if (!location.protocol.startsWith('http')) { setStatus(fileNote, true); return; }
	const btn = $('run');
	btn.disabled = true;
	btn.textContent = 'Running…';
	setStatus('<span class="spin"></span>Running the simulation…');
	const t0 = performance.now();
	try {
		const res = await fetch('/api/run?' + query());
		const body = await res.json().catch(() => ({ error: `server answered ${res.status}` }));
		if (!res.ok || body.error) throw new Error(body.error || `server answered ${res.status}`);
		show(body);
		const secs = (performance.now() - t0) / 1000;
		setStatus(`Done in ${secs.toFixed(2)} s, including the classical comparisons and sending ${data.layers.length} layers to the page.`);
	} catch (e) {
		setStatus('Error: ' + e.message, true);
	} finally {
		btn.disabled = false;
		btn.textContent = 'Run';
	}
}

function fillControls(p) {
	F('problem').value = p.problem;
	syncControls();
	for (const k of ['n', 'layers', 'dt', 'tilt', 'seed']) F(k).value = p[k];
	for (const k of Object.keys(EXTRAS)) if (p[k] !== undefined && p[k] !== null) F(k).value = p[k];
	for (const k of ['engine', 'schedule', 'init']) F(k).value = p[k];
	syncControls();
}

// show the inputs of the chosen problem; keep the controls to combinations the simulator accepts
function syncControls() {
	const problem = F('problem').value, rank1 = problem === 'partition';
	for (const l of form.querySelectorAll('label[data-problem]')) l.hidden = l.dataset.problem !== problem;
	F('engine').querySelector('[value=hw]').disabled = !rank1;
	if (!rank1) F('engine').value = 'float';
	F('schedule').querySelector('[value=lut]').disabled = !rank1;
	if (!rank1) F('schedule').value = 'ramp';
	const lut = F('schedule').value === 'lut';
	F('layers').max = lut ? 4 : 5000;
	if (lut && +F('layers').value > 4) F('layers').value = 4;
	F('dt').disabled = lut;
	F('n').max = MAX_N[problem] || 12000;
}


// ---------------------------------------------------------------
// studies: form, start, poll, cancel
// ---------------------------------------------------------------

function syncStudyForm() {
	const st = SF('study').value;
	for (const o of SF('problem').querySelectorAll('[data-sweep-only]')) o.disabled = st === 'scaling';
	if (st === 'scaling' && SF('problem').selectedOptions[0].disabled) SF('problem').value = 'sk';
	const problem = SF('problem').value;
	for (const l of studyForm.querySelectorAll('label[data-study]')) {
		const inStudy = l.dataset.study.split(' ').includes(st);
		l.hidden = !inStudy || (l.dataset.problem !== undefined && l.dataset.problem !== problem);
	}
}

function applyPreset(which) {
	const st = SF('study').value, p = STUDY_PRESETS[st][which];
	for (const [k, v] of Object.entries(p)) if (SF(k)) SF(k).value = v;
	syncStudyForm();
	for (const b of $('presets').querySelectorAll('button')) b.classList.toggle('on', b.dataset.preset === which);
}

function studyQuery() {
	const parts = [`study=${SF('study').value}`];
	for (const l of studyForm.querySelectorAll('label')) {
		if (l.hidden || l.classList.contains('file')) continue;
		const el = l.querySelector('input, select');
		if (!el || !el.name || el.name === 'study' || el.value === '') continue;
		parts.push(`${el.name}=${encodeURIComponent(el.value.replace(/\s+/g, ''))}`);
	}
	return parts.join('&');
}

let job = null;		// { id, t0, study, timer }

function studyMsg(html, error = false) {
	$('studyMessage').innerHTML = html;
	$('studyMessage').classList.toggle('error', error);
}

function setProgress(p) {
	$('studyFill').style.width = (100 * Math.max(0, Math.min(1, p))).toFixed(1) + '%';
}

function setRunningUi(on) {
	studyRunning = on;
	$('studyStart').disabled = on;
	$('studyStart').textContent = on ? 'Running…' : 'Start';
	$('studyCancel').disabled = !on;
}

async function startStudy() {
	if (!location.protocol.startsWith('http')) { studyMsg(fileNote, true); return; }
	if (job) return;
	setRunningUi(true);
	setProgress(0);
	studyMsg('<span class="spin"></span>Starting…');
	try {
		const res = await fetch('/api/study/start?' + studyQuery());
		const body = await res.json().catch(() => ({ error: `server answered ${res.status}` }));
		if (!res.ok || body.error) throw new Error(body.error || `server answered ${res.status}`);
		job = { id: body.id, t0: performance.now(), study: SF('study').value };
		poll(job);
	} catch (e) {
		studyMsg('Error: ' + e.message, true);
		setRunningUi(false);
	}
}

async function poll(j) {
	if (job !== j) return;
	try {
		const res = await fetch(`/api/study/status?id=${j.id}`);
		const s = await res.json().catch(() => ({ error: `server answered ${res.status}` }));
		if (job !== j) return;
		if (!res.ok) throw new Error(s.error || `server answered ${res.status}`);
		const elapsed = (performance.now() - j.t0) / 1000;
		setProgress(s.progress);
		if (s.done) {
			job = null;
			setRunningUi(false);
			if (s.result) {
				setProgress(1);
				showStudy(s.result);
				studyMsg(`Done in ${fmtTime(elapsed)}.`);
			} else {
				studyMsg(s.error === 'cancelled' ? 'Cancelled.' : 'Error: ' + s.error, s.error !== 'cancelled');
			}
			return;
		}
		// progress is by estimated work, so the time left is roughly proportional
		const p = s.progress;
		const eta = p > 0.02 && elapsed > 1 ? ` · about ${fmtTime(elapsed * (1 - p) / p)} left` : '';
		studyMsg(`<span class="spin"></span>${(100 * p).toFixed(0)}% · ${s.message || 'running'} · ${fmtTime(elapsed)} elapsed${eta}`);
		if (s.partial) showStudy(s.partial, true);
		j.timer = setTimeout(() => poll(j), 500);
	} catch (e) {
		if (job === j) job = null;
		setRunningUi(false);
		studyMsg('Error: ' + e.message, true);
	}
}

async function cancelStudy() {
	if (!job) return;
	studyMsg('<span class="spin"></span>Cancelling…');
	try { await fetch(`/api/study/cancel?id=${job.id}`); } catch (e) { /* the next poll reports it */ }
}

function fillStudyForm(d) {
	SF('study').value = d.study;
	const p = d.params || {};
	if (d.study === 'scaling') {
		SF('problem').value = p.problem || 'sk';
		SF('sizes').value = (d.sizes || []).join(',');
		SF('seeds').value = d.seeds;
		SF('layers').value = p.layers;
		SF('dt').value = p.dt;
	} else if (d.study === 'sweep') {
		SF('problem').value = p.problem || 'sk';
		SF('n').value = p.n;
		for (const k of Object.keys(EXTRAS)) if (p[k] !== undefined && p[k] !== null) SF(k).value = p[k];
		SF('layers_list').value = d.layers_list.join(',');
		SF('dt_list').value = d.dt_list.join(',');
		SF('seeds').value = d.seeds;
	} else {
		SF('n').value = d.n;
		SF('trajectories').value = d.trajectories;
		SF('times').value = d.times.length;
		SF('t_max').value = d.times[d.times.length - 1];
	}
	if (p.seed !== undefined) SF('seed').value = p.seed;
	syncStudyForm();
	for (const b of $('presets').querySelectorAll('button')) b.classList.remove('on');
}


// ---------------------------------------------------------------
// studies: drawing
// ---------------------------------------------------------------

function showStudy(d, running = false) {
	study = d;
	study._running = running;
	if (mode === 'studies') renderStudy();
}

function renderStudy() {
	const d = study;
	$('sweepView').hidden = !d || d.study !== 'sweep';
	$('studyCard2').hidden = !d || d.study !== 'squeezing';
	$('studies').classList.toggle('single', $('studyCard2').hidden);
	if (!d) {
		$('studyAboutTitle').textContent = 'About the studies';
		$('studyAbout').innerHTML = studiesIntro();
		$('studyChartTitle').textContent = 'Result';
		$('studyChartCaption').textContent = 'Start a study, or load a study\'s JSON file, to see its result here.';
		const { ctx, w, h } = setup($('studyChart'));
		ctx.fillStyle = T.muted;
		ctx.textAlign = 'center';
		ctx.textBaseline = 'middle';
		ctx.fillText('No study yet', w / 2, h / 2);
		$('studyChart')._geom = null;
		$('studyLegend').innerHTML = '';
		return;
	}
	if (d.study === 'scaling') renderScaling(d);
	else if (d.study === 'sweep') renderSweep(d);
	else renderSqueezing(d);
}

function studiesIntro() {
	return `<p>A study runs many simulations and draws the results as curves. It runs in the background on the server; the charts fill in as results arrive, and Cancel stops it.</p>
		<ul>
			<li><b>Spin glass toward N → ∞:</b> SK instances at growing N, extrapolated to infinite size and compared with the Parisi ground-state energy.</li>
			<li><b>Schedule sweep:</b> the final result over a grid of layer counts and ramp sizes dt, the tuning of γ and β that mean-field QAOA needs.</li>
			<li><b>Squeezing from classical runs:</b> one-axis twisting, where averaging many classical runs recovers a quantum effect that one mean-field run misses.</li>
		</ul>
		<p><b>Quick</b> takes seconds; <b>Big</b> about a minute on a laptop.</p>`;
}

function timingTiles(d, note) {
	const ratio = d.cpu_seconds / d.fpga_seconds;
	return `<div class="subhead">Time</div><div class="tiles">` +
		tile(`CPU: ${d.cpu}`, fmtTime(d.cpu_seconds), `simulation time; ${fmtTime(d.wall_seconds)} wall clock`) +
		tile('FPGA estimate', fmtTime(d.fpga_seconds), 'the same runs') +
		tile(ratio >= 1 ? 'FPGA faster by' : 'CPU faster by', `${fmt(ratio >= 1 ? ratio : 1 / ratio, 3)}×`, 'simulation time') +
		`</div><p class="note">${note}</p>` +
		(d.summary ? `<div class="subhead">Summary</div><pre class="summary">${d.summary.replace(/&/g, '&amp;').replace(/</g, '&lt;')}</pre>` : '');
}

function renderScaling(d) {
	const pts = d.points || [];
	const prob = d.params && d.params.problem === 'sk-gauss' ? 'Gaussian' : '±1/√N';
	$('studyAboutTitle').textContent = 'Spin glass toward N → ∞';
	$('studyAbout').innerHTML = `<p>Mean-field QAOA on Sherrington–Kirkpatrick spin glasses (couplings ${prob}) at sizes N = ${(d.sizes || []).map((n) => n.toLocaleString()).join(', ')}, ${d.seeds} random instances each, ${d.params.layers} layers with dt = ${d.params.dt}.</p>
		<p>The energy per spin of a finite system differs from the infinite one by roughly a constant times N<sup>−2/3</sup>. Plotted against x = N<sup>−2/3</sup>, the points fall on a line, and where the line meets x = 0 is the estimate for N → ∞ (the left edge of the chart).</p>
		<p>The dashed line at −0.7632 is the exact ground-state energy per spin for N → ∞ (Parisi). The gap between it and the extrapolated mean-field value is what this schedule can't reach, however large the system; more layers narrow it.</p>
		<p>A dense J is the proposed Agilex 7 design's job, so this is where an FPGA would pay off: one N = 8,000 layer is about 0.6 ms on a laptop and an estimated 20 µs on the FPGA.</p>`;
	$('studyChartTitle').textContent = 'Energy per spin against N^(−2/3)';
	$('studyChartCaption').innerHTML = `Each point is the mean over the instances at one N, with its standard error; faint dots are single instances. Dashed lines are fits, extended to N → ∞ at x = 0.${d._running ? ' <b>Still running</b>: more sizes are on the way.' : ''}`;
	const xs = pts.map((p) => p.x), xmax = pts.length ? Math.max(...xs) * 1.08 : 0.03;
	const se = (s, c) => s / Math.sqrt(Math.max(1, c));
	const series = [
		{ name: 'single instances', color: T.a, xs: pts.flatMap((p) => p.values.map(() => p.x)), ys: pts.flatMap((p) => p.values), dots: true, r: 2.2, alpha: 0.3 },
		{ name: 'mean field, mean ± s.e.', color: T.a, xs, ys: pts.map((p) => p.mean), err: pts.map((p) => se(p.std, p.count)), dots: true, r: 4.5 },
		{ name: '+ local search', color: T.c3, xs, ys: pts.map((p) => p.local_mean), err: pts.map((p) => se(p.local_std, p.count)), dots: true, r: 4.5 },
	];
	const pm = (f) => `${f.e_inf.toFixed(4)} ± ${f.e_inf_err.toFixed(4)}`;
	const fits = [];
	// the lower of the two extrapolations gets its label below, the other above
	const mfLower = d.fit && d.fit_local && d.fit.e_inf < d.fit_local.e_inf;
	if (d.fit) fits.push({ color: T.a, a: d.fit.e_inf, b: d.fit.slope, x0: 0, x1: xmax, label: pm(d.fit), below: mfLower });
	if (d.fit_local) fits.push({ color: T.c3, a: d.fit_local.e_inf, b: d.fit_local.slope, x0: 0, x1: xmax, label: pm(d.fit_local), below: !mfLower });
	const draw = () => xyChart($('studyChart'), {
		series, fits, xmin: 0, xmax,
		hlines: [{ label: `ground state for N → ∞ (Parisi), ${d.parisi}`, value: d.parisi }],
		topTicks: pts.map((p) => ({ x: p.x, label: 'N=' + p.n.toLocaleString() })),
		xTitle: 'x = N^(−2/3)    (N → ∞ at x = 0)',
		xFmt: (v) => fmt(v, 2),
		pickXs: xs,
		tip: (i) => {
			const p = pts[i];
			return `<b>N = ${p.n.toLocaleString()}</b><div>${p.count} instances</div><div><span class="sw dot" style="background:${T.a}"></span>mean field ${fmt(p.mean, 4)} ± ${fmt(p.std, 2)}</div><div><span class="sw dot" style="background:${T.c3}"></span>+ local search ${fmt(p.local_mean, 4)} ± ${fmt(p.local_std, 2)}</div>`;
		},
		render: draw,
	});
	draw();
	$('studyLegend').innerHTML = key(T.a, 'mean field (rounded), mean ± standard error', 'dot') + key(T.c3, 'mean field + local search', 'dot') +
		key(T.a, 'single instances', 'dot faint') + key(null, 'fits and the Parisi energy', 'dash');

	let html = '';
	if (d.fit) {
		html += '<div class="tiles">' +
			tile('Mean field, N → ∞', fmt(d.fit.e_inf, 4), `± ${d.fit.e_inf_err.toFixed(4)}`) +
			(d.fit_local ? tile('+ local search, N → ∞', fmt(d.fit_local.e_inf, 4), `± ${d.fit_local.e_inf_err.toFixed(4)}`) : '') +
			tile('Ground state, N → ∞', fmt(d.parisi, 4), 'Parisi') +
			tile('Gap to the ground state', fmt(d.fit.e_inf - d.parisi, 3), 'per spin, mean field') +
			'</div>';
	}
	html += `<table><thead><tr><th class="num">N</th><th class="num">Instances</th><th class="num">Mean field</th><th class="num">+ local search</th><th class="num">CPU time</th></tr></thead><tbody>` +
		pts.map((p) => `<tr><td class="num">${p.n.toLocaleString()}</td><td class="num">${p.count}</td><td class="num">${fmt(p.mean, 4)} ± ${fmt(p.std, 2)}</td><td class="num">${fmt(p.local_mean, 4)} ± ${fmt(p.local_std, 2)}</td><td class="num">${fmtTime(p.seconds)}</td></tr>`).join('') +
		'</tbody></table>';
	html += timingTiles(d, 'FPGA estimate: the dense-J Agilex 7 design (not built), J in block RAM up to N = 8,700 and in HBM above.');
	$('studyResults').innerHTML = html;
}

function renderSweep(d) {
	const view = sweepView;
	for (const b of $('sweepView').querySelectorAll('button')) b.classList.toggle('on', b.dataset.view === view);
	const grid = d[view], L = d.layers_list, D = d.dt_list;
	const energy = view === 'grid_energy';
	const higher = energy ? false : d.metric.higher_is_better;
	const fraction = !energy && /^Fraction/.test(d.metric.title);
	const title = energy ? 'Energy per spin' : d.metric.title + (view === 'grid_local' ? ', + local search' : ', mean field');
	const prob = d.params.problem;
	$('studyAboutTitle').textContent = 'Schedule sweep';
	$('studyAbout').innerHTML = `<p>${probName(prob)} with N = ${d.params.n.toLocaleString()}, ${d.seeds} random instances per setting. Each cell is the mean result for one number of layers and one ramp size dt (γ rises from 0 to dt and β falls from dt to 0).</p>
		<p>Mean-field QAOA needs its schedule tuned: too few layers or too large a dt and the spins can't follow the changing field; too small a dt and the ramp never gets far enough. Slower ramps (more layers) generally do better, and the best dt shifts as layers change.</p>
		<p>This is also why an FPGA with several pipelines would help: each cell is independent, so several (γ, β) settings could run side by side.</p>`;
	$('studyChartTitle').textContent = title;
	$('studyChartCaption').innerHTML = `Rows: layers; columns: dt. Stronger blue is better (${higher ? 'higher' : 'lower'} ${energy ? 'energy' : d.metric.title.toLowerCase()}): cells are shaded by rank, so the best region stands out even when the values are close. The best cell is outlined.${d._running ? ' <b>Still running</b>: blank cells are on the way.' : ''}`;

	// shade by rank among the finished cells, worst 0 to best 1
	const vals = [];
	for (const row of grid) for (const v of row) if (Number.isFinite(v)) vals.push(v);
	vals.sort((a, b) => (higher ? a - b : b - a));
	const good = (v) => (vals.length > 1 ? vals.indexOf(v) / (vals.length - 1) : 1);
	let best = null;
	grid.forEach((row, i) => row.forEach((v, j) => {
		if (v === null || !isFinite(v)) return;
		if (!best || (higher ? v > best.v : v < best.v)) best = { i, j, v };
	}));

	const cv = $('studyChart');
	const draw = () => {
		cv.style.height = Math.max(240, Math.min(420, 44 + L.length * 40)) + 'px';
		const { ctx, w, h } = setup(cv);
		const pad = { l: 58, r: 10, t: 8, b: 40 };
		const pw = w - pad.l - pad.r, ph = h - pad.t - pad.b;
		const cw = pw / D.length, chh = ph / L.length;
		const hv = hover[cv.id];
		const cellY = (i) => pad.t + (L.length - 1 - i) * chh;	// fewest layers at the bottom
		ctx.font = '11px ' + getComputedStyle(document.body).fontFamily;
		grid.forEach((row, i) => row.forEach((v, j) => {
			const x = pad.l + j * cw, y = cellY(i);
			if (v === null || !isFinite(v)) {
				ctx.save();
				ctx.setLineDash([3, 3]);
				ctx.strokeStyle = T.axis;
				ctx.strokeRect(x + 2.5, y + 2.5, cw - 5, chh - 5);
				ctx.restore();
				return;
			}
			const t = good(v);
			ctx.fillStyle = css(mixRgb(T.mid, T.a, 0.06 + 0.94 * t));
			ctx.fillRect(x + 1, y + 1, cw - 2, chh - 2);
			if (cw >= 46 && chh >= 18) {
				ctx.textAlign = 'center';
				ctx.textBaseline = 'middle';
				ctx.fillStyle = t > 0.55 ? '#ffffff' : T.text;
				ctx.fillText(fraction ? pct(v, 1) : fmt(v, 3), x + cw / 2, y + chh / 2);
			}
		}));
		if (best) {
			ctx.strokeStyle = T.text;
			ctx.lineWidth = 2.5;
			ctx.strokeRect(pad.l + best.j * cw + 1.5, cellY(best.i) + 1.5, cw - 3, chh - 3);
		}
		if (hv !== null && hv !== undefined) {
			ctx.strokeStyle = T.text2;
			ctx.lineWidth = 1.5;
			ctx.strokeRect(pad.l + (hv % D.length) * cw + 0.5, cellY(Math.floor(hv / D.length)) + 0.5, cw - 1, chh - 1);
		}
		ctx.fillStyle = T.muted;
		ctx.textAlign = 'right';
		ctx.textBaseline = 'middle';
		L.forEach((l, i) => ctx.fillText(String(l), pad.l - 6, cellY(i) + chh / 2));
		ctx.textAlign = 'center';
		ctx.textBaseline = 'top';
		D.forEach((v, j) => ctx.fillText(String(v), pad.l + (j + 0.5) * cw, pad.t + ph + 6));
		ctx.fillStyle = T.text2;
		ctx.fillText('dt (ramp size)', pad.l + pw / 2, pad.t + ph + 22);
		ctx.save();
		ctx.translate(12, pad.t + ph / 2);
		ctx.rotate(-Math.PI / 2);
		ctx.fillText('layers', 0, 0);
		ctx.restore();
		cv._geom = {
			render: draw,
			pick: (px, py) => {
				const j = Math.floor((px - pad.l) / cw), r = Math.floor((py - pad.t) / chh);
				const i = L.length - 1 - r;
				return j >= 0 && j < D.length && i >= 0 && i < L.length ? i * D.length + j : null;
			},
			tip: (k) => {
				const i = Math.floor(k / D.length), j = k % D.length;
				const g = (name) => d[name][i][j];
				if (g('grid') === null) return `<b>${L[i]} layers, dt ${D[j]}</b><div>still running</div>`;
				return `<b>${L[i]} layers, dt ${D[j]}</b><div>${d.metric.title}: ${fmt(g('grid'), 4)}</div><div>+ local search: ${fmt(g('grid_local'), 4)}</div><div>energy per spin: ${fmt(g('grid_energy'), 4)}</div>`;
			},
		};
	};
	draw();
	const c0 = css(mixRgb(T.mid, T.a, 0.06)), c1 = css(mixRgb(T.mid, T.a, 1));
	$('studyLegend').innerHTML = vals.length
		? `<span class="key">worst ${fmt(vals[0], 3)}<span class="ramp" style="background:linear-gradient(90deg, ${c0}, ${c1})"></span>${fmt(vals[vals.length - 1], 3)} best</span><span>outlined: best</span><span>dashed: not run yet</span>`
		: '<span>dashed: not run yet</span>';

	let html = '';
	if (best) {
		html += '<div class="tiles">' +
			tile('Best setting', `${L[best.i]} layers, dt ${D[best.j]}`, title) +
			tile('Its value', fmt(best.v, 4), `mean over ${d.seeds} instances`) +
			'</div>';
	}
	html += timingTiles(d, 'FPGA estimate: the same runs on the FPGA designs described for single runs (most are estimates for designs that aren\'t built).');
	$('studyResults').innerHTML = html;
}

function probName(p) {
	return {
		partition: 'Number partitioning', sk: 'The SK spin glass (±1/√N)', 'sk-gauss': 'The SK spin glass (Gaussian)',
		maxcut: 'Max-Cut on a random regular graph', chain: 'The long-range chain', hopfield: 'The Hopfield memory',
	}[p] || p;
}

function renderSqueezing(d) {
	const ts = d.times, n = d.n;
	$('studyAboutTitle').textContent = 'Squeezing from classical runs';
	$('studyAbout').innerHTML = `<p><b>One-axis twisting.</b> ${n.toLocaleString()} spins all start along +x and evolve under H = χS<sub>z</sub>², the all-to-all coupling with every weight equal: exactly what <code>VIRAL_n_spin</code> computes, with no mixer. The quantum state becomes squeezed: the spread of the total spin perpendicular to its mean drops below the coherent-state value. The squeezing parameter ξ² = N·(smallest perpendicular variance)/⟨S<sub>x</sub>⟩² measures it; below 0 dB is squeezed.</p>
		<p><b>Mean field</b> (a single classical run) gives no squeezing at all: from exactly +x every field is zero and nothing moves.</p>
		<p><b>Many classical runs</b> recover it: each of the ${d.trajectories.toLocaleString()} runs starts every spin at (1, ±1, ±1) with random signs (the discrete truncated Wigner method), and averaging the runs gives the quantum fluctuations. The exact curve comes from the symmetric (Dicke) states, which is possible here because all spins are equivalent.</p>`;
	$('studyChartTitle').textContent = 'Squeezing ξ² against χt';
	$('studyChartCaption').innerHTML = `In dB: 10 log<sub>10</sub> ξ², lower is more squeezed. Rings mark the best squeezing.${d._running ? ` <b>Still running</b>: ${d.done.toLocaleString()} of ${d.trajectories.toLocaleString()} trajectories so far.` : ''}`;
	const ex = d.exact.xi2_db, dt = d.dtwa.xi2_db;
	// keep the y axis on the interesting part: from the best squeezing to a bit above 0 dB
	const deepest = Math.min(0, ...ex.filter(Number.isFinite));
	const cap = Math.max(15, -deepest * 0.6);
	const lim = (v) => (Number.isFinite(v) && v <= cap ? v : null);
	const draw = () => xyChart($('studyChart'), {
		series: [
			{ name: 'exact', color: T.a, xs: ts, ys: ex.map(lim), line: true },
			{ name: 'classical runs (DTWA)', color: T.b, xs: ts, ys: dt.map(lim), dots: true, r: 2.6 },
		],
		hlines: [{ label: 'mean field: one classical run, no squeezing (0 dB)', value: 0 }],
		marks: [
			{ x: d.best_exact.t, y: d.best_exact.xi2_db, color: T.a, label: `exact ${fmt(d.best_exact.xi2_db, 3)} dB` },
			{ x: d.best_dtwa.t, y: d.best_dtwa.xi2_db, color: T.b, label: `DTWA ${fmt(d.best_dtwa.xi2_db, 3)} dB`, above: true },
		],
		xTitle: 'χt', xFmt: (v) => fmt(v, 2), yFmt: (v) => fmt(v, 3) + ' dB',
		pickXs: ts,
		tip: (i) => `<b>χt = ${fmt(ts[i], 4)}</b><div><span class="sw line" style="background:${T.a}"></span>exact ${fmt(ex[i], 4)} dB</div><div><span class="sw dot" style="background:${T.b}"></span>DTWA ${fmt(dt[i], 4)} dB</div>`,
		render: draw,
	});
	draw();
	$('studyLegend').innerHTML = key(T.a, 'exact (Dicke states)', 'line') + key(T.b, `classical runs (DTWA), ${d.done.toLocaleString()} averaged`, 'dot') + key(null, 'mean field', 'dash');

	$('studyChart2Title').textContent = 'Mean spin ⟨S_x⟩ / (N/2)';
	$('studyChart2Caption').textContent = 'The mean spin shrinks as the state spreads around the equator; past the best squeezing it falls toward 0 and ξ² grows.';
	const exs = d.exact.sx, dts = d.dtwa.sx;
	const draw2 = () => xyChart($('studyChart2'), {
		series: [
			{ name: 'exact', color: T.a, xs: ts, ys: exs, line: true },
			{ name: 'DTWA', color: T.b, xs: ts, ys: dts, dots: true, r: 2.6 },
		],
		hlines: [{ label: 'mean field: stays at 1', value: 1 }],
		xTitle: 'χt', xFmt: (v) => fmt(v, 2),
		pickXs: ts,
		tip: (i) => `<b>χt = ${fmt(ts[i], 4)}</b><div><span class="sw line" style="background:${T.a}"></span>exact ${fmt(exs[i], 4)}</div><div><span class="sw dot" style="background:${T.b}"></span>DTWA ${fmt(dts[i], 4)}</div>`,
		render: draw2,
	});
	draw2();
	$('studyLegend2').innerHTML = key(T.a, 'exact', 'line') + key(T.b, 'classical runs (DTWA)', 'dot') + key(null, 'mean field', 'dash');

	let html = '<div class="tiles">' +
		tile('Best squeezing, exact', `${fmt(d.best_exact.xi2_db, 3)} dB`, `at χt = ${fmt(d.best_exact.t, 3)}`) +
		tile('Best squeezing, DTWA', `${fmt(d.best_dtwa.xi2_db, 3)} dB`, `at χt = ${fmt(d.best_dtwa.t, 3)}`) +
		tile('Mean field', '0 dB', 'no squeezing') +
		tile('Trajectories', d.done.toLocaleString(), `of ${d.trajectories.toLocaleString()}`) +
		'</div>';
	html += timingTiles(d, 'FPGA estimate: one VIRAL_n_spin-style pipeline (one spin per clock at 300 MHz) running every trajectory and time point. This is rank-1 work, so 8 CPU cores with SIMD keep up with it.');
	$('studyResults').innerHTML = html;
}


// ---------------------------------------------------------------
// mode, files, wiring
// ---------------------------------------------------------------

function setMode(m) {
	mode = m;
	const studies = m === 'studies';
	for (const b of document.querySelectorAll('.tabs button')) b.setAttribute('aria-selected', String(b.dataset.tab === m));
	form.hidden = studies;
	$('status').hidden = studies;
	$('layerbar').hidden = studies;
	$('grid').hidden = studies;
	studyForm.hidden = !studies;
	$('studyStatus').hidden = !studies;
	$('studies').hidden = !studies;
	if (studies) setPlaying(false);
	hideTip();
	renderAll();
}

async function loadFile(e) {
	const file = e.target.files[0];
	if (!file) return;
	const msg = (html, err) => (mode === 'studies' ? studyMsg(html, err) : setStatus(html, err));
	try {
		let d = JSON.parse(await file.text());
		if (d.result && d.result.study) d = d.result;		// a saved /api/study/status answer
		if (d.study) {
			setMode('studies');
			fillStudyForm(d);
			showStudy(d);
			studyMsg(`Loaded ${file.name}.`);
		} else if (d.layers && d.problem) {
			setMode('run');
			show(d, true);
			setStatus(`Loaded ${file.name}.`);
		} else throw new Error('not a viral-sim run or study');
	} catch (err) {
		msg('Could not load the file: ' + err.message, true);
	}
	e.target.value = '';
}

form.addEventListener('submit', (e) => { e.preventDefault(); run(); });

F('problem').addEventListener('change', () => {
	const d = DEFAULTS[F('problem').value];
	F('n').value = d.n;
	F('layers').value = d.layers;
	F('dt').value = d.dt;
	F('schedule').value = 'ramp';
	F('engine').value = F('problem').value === 'partition' ? 'hw' : 'float';
	syncControls();
});
F('schedule').addEventListener('change', syncControls);

studyForm.addEventListener('submit', (e) => { e.preventDefault(); startStudy(); });
SF('study').addEventListener('change', () => applyPreset('quick'));
SF('problem').addEventListener('change', () => {
	if (SF('study').value === 'sweep') SF('n').value = DEFAULTS[SF('problem').value].n;
	syncStudyForm();
	for (const b of $('presets').querySelectorAll('button')) b.classList.remove('on');
});
$('presets').addEventListener('click', (e) => {
	const b = e.target.closest('button');
	if (b) applyPreset(b.dataset.preset);
});
$('studyCancel').addEventListener('click', cancelStudy);
$('sweepView').addEventListener('click', (e) => {
	const b = e.target.closest('button');
	if (!b || !study) return;
	sweepView = b.dataset.view;
	renderStudy();
});

for (const inp of document.querySelectorAll('.loadfile')) inp.addEventListener('change', loadFile);
for (const b of document.querySelectorAll('.tabs button')) b.addEventListener('click', () => setMode(b.dataset.tab));

$('layer').addEventListener('input', (e) => { setPlaying(false); setLayer(+e.target.value); });
$('play').addEventListener('click', () => setPlaying(!playing));

$('setsMode').addEventListener('click', (e) => {
	const b = e.target.closest('button');
	if (!b || b.disabled || !data) return;
	setsMode = b.dataset.mode;
	renderProblem();
});
$('skView').addEventListener('click', (e) => {
	const b = e.target.closest('button');
	if (!b || !data) return;
	skView = b.dataset.view;
	renderProblem();
});

for (const id of ['problemCanvas', 'problemCanvas2', 'energy', 'metric', 'schedule', 'hist', 'studyChart', 'studyChart2']) attachHover($(id));
attachSphere();

let resizeQueued = false;
const ro = new ResizeObserver(() => {
	if (resizeQueued) return;
	resizeQueued = true;
	requestAnimationFrame(() => { resizeQueued = false; renderAll(); });
});
ro.observe($('grid'));
ro.observe($('studies'));

matchMedia('(prefers-color-scheme: dark)').addEventListener('change', renderAll);

syncControls();
applyPreset('quick');
renderAll();
if (location.protocol.startsWith('http')) run();
else setStatus('Opened as a file: load a JSON file from <code>viral-sim run --json</code> or <code>viral-sim study --json</code>, or use <code>viral-sim serve</code> to run simulations.');
