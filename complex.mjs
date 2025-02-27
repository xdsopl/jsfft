/*
Complex math

Copyright 2025 Ahmet Inan <xdsopl@gmail.com>
*/

export default class Complex {
	constructor(a, b) {
		this.real = a || 0;
		this.imag = b || 0;
	}

	set(a, b) {
		this.real = a;
		this.imag = b;
		return this;
	}

	copy(other) {
		this.real = other.real;
		this.imag = other.imag;
		return this;
	}

	norm() {
		return this.real * this.real + this.imag * this.imag;
	}

	abs() {
		return Math.sqrt(this.norm());
	}

	arg() {
		return Math.atan2(this.imag, this.real);
	}

	polar(a, b) {
		this.real = a * Math.cos(b);
		this.imag = a * Math.sin(b);
		return this;
	}

	conj() {
		this.imag = -this.imag;
		return this;
	}

	add(other) {
		this.real += other.real;
		this.imag += other.imag;
		return this;
	}

	sub(other) {
		this.real -= other.real;
		this.imag -= other.imag;
		return this;
	}

	scale(s) {
		this.real *= s;
		this.imag *= s;
		return this;
	}

	mul(other) {
		const tmp = this.real * other.real - this.imag * other.imag;
		this.imag = this.real * other.imag + this.imag * other.real;
		this.real = tmp;
		return this;
	}

	divide(s) {
		this.real /= s;
		this.imag /= s;
		return this;
	}

	div(other) {
		const den = other.norm();
		const tmp = (this.real * other.real + this.imag * other.imag) / den;
		this.imag = (this.imag * other.real - this.real * other.imag) / den;
		this.real = tmp;
		return this;
	}

	toString() {
		return `(${this.real} + ${this.imag}i)`;
	}
}
