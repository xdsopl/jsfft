/*
Some helper functions

Copyright 2025 Ahmet Inan <xdsopl@gmail.com>
*/

export function hann(n, N) {
	return 0.5 * (1.0 - Math.cos((2.0 * Math.PI * n) / (N - 1)));
}

export function sinc(x) {
	if (x === 0)
		return 1;
	x *= Math.PI;
	return Math.sin(x) / x;
}

export function lowPass(cutoff, rate, n, N) {
	const f = 2 * cutoff / rate;
	const x = n - (N - 1) / 2;
	return f * sinc(f * x);
}

