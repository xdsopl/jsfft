/*
Frequency sweep test of STFT

Copyright 2025 Ahmet Inan <xdsopl@gmail.com>
*/

import Complex from './complex.mjs';
import ShortTimeFourierTransform from './stft.mjs';

const length = 420;
const stft = new ShortTimeFourierTransform(length, 3);
const input = new Complex();
const factor = 10000, range = length * factor;
function dB(power) {
	return 10.0 * Math.log10(power);
}
for (let acc = 0, freq = -10 * factor; freq <= 10 * factor; acc = (acc + ++freq) % range)
	if (stft.push(input.polar(1, (acc * 2 * Math.PI) / range)))
		console.log(freq / factor + " "
			+ dB(stft.power[length - 1]) + " " + dB(stft.power[0]) + " " + dB(stft.power[1]));

