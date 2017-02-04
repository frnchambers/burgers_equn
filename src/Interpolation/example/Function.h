

double f ( double x ) {
	return (x*x*x*x) - x*x + 1;
//	return x*x*x*x*x;
//	return std::exp(x);
//	return std::cos(x);
//	return 3.0*x;
}

double df ( double x ) {
	return 4.0*x*x*x - 2.0*x;
//	return 5.0*x*x*x*x;
//	return std::exp(x);
//	return -std::sin(x);
//	return 3.0;
}

double d2f ( double x ) {
	return 12.0*x*x - 2.0;
//	return 20.0*x*x*x;
//	return std::exp(x);
//	return -std::cos(x);
//	return 0.0;
}
