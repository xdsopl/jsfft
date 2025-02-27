/*
Short Time Fourier Transform

Copyright 2025 Ahmet Inan <xdsopl@gmail.com>
*/

import FastFourierTransform from './fft.mjs';
import Complex from './complex.mjs';
import { hann, lowPass } from './utils.mjs';

export default class ShortTimeFourierTransform {
	constructor(length, overlap) {
		this.fft = new FastFourierTransform(length);
		this.prev = Array(length * overlap);
		for (let i = 0; i < length * overlap; ++i)
			this.prev[i] = new Complex();
		this.fold = new Array(length);
		for (let i = 0; i < length; ++i)
			this.fold[i] = new Complex();
		this.freq = new Array(length);
		for (let i = 0; i < length; ++i)
			this.freq[i] = new Complex();
		this.temp = new Complex();
		this.power = new Array(length);
		this.weight = new Array(length * overlap);
		for (let i = 0; i < length * overlap; ++i)
			this.weight[i] = lowPass(1, length, i, length * overlap) * hann(i, length * overlap);
		this.index = 0;
	}

	push(input) {
		this.prev[this.index].copy(input);
		this.index = (this.index + 1) % this.prev.length;
		if (this.index % this.fold.length != 0)
			return false;
		for (let i = 0; i < this.fold.length; ++i, this.index = (this.index + 1) % this.prev.length)
			this.fold[i].copy(this.prev[this.index]).scale(this.weight[i]);
		for (let i = this.fold.length; i < this.prev.length; ++i, this.index = (this.index + 1) % this.prev.length)
			this.fold[i % this.fold.length].add(this.temp.copy(this.prev[this.index]).scale(this.weight[i]));
		this.fft.forward(this.freq, this.fold);
		for (let i = 0; i < this.power.length; ++i)
			this.power[i] = this.freq[i].norm();
		return true;
	}
}