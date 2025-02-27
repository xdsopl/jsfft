/*
Test bench

Copyright 2025 Ahmet Inan <xdsopl@gmail.com>
*/

import FastFourierTransform from './fft.mjs';
import Complex from './complex.mjs';

const length = parseInt(process.argv[2], 10);
const fft = new FastFourierTransform(length);
const buf0 = new Array(length);
const buf1 = new Array(length);
const buf2 = new Array(length);
for (let i = 0; i < length; ++i) {
	const a = Math.random() * 2 - 1;
	const b = Math.random() * 2 - 1;
	buf0[i] = new Complex(a, b);
}
for (let i = 0; i < length; ++i)
	buf1[i] = new Complex();
for (let i = 0; i < length; ++i)
	buf2[i] = new Complex();
fft.forward(buf1, buf0);
fft.backward(buf2, buf1);
const factor = 1.0 / length;
for (let i = 0; i < length; ++i)
	buf2[i].scale(factor);
let maxError = 0;
for (let i = 0; i < length; ++i)
	maxError = Math.max(maxError, buf2[i].sub(buf0[i]).abs());
console.log(`max error = ${maxError}`);
if (maxError > 1.0e-14)
	throw new Error(`Transform error too high: ${maxError}`);
