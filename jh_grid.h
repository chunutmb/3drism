/*This file is for code to manipulate a grid vector or what ever*/

#ifndef JH_SHIFTING
  #define JH_SHIFTING

	void   shift_origin_complex( fftw_complex *, fftw_complex *, ENV_PAR  );	/*in, out, var*/
	void unshift_origin_complex( fftw_complex *, fftw_complex *, ENV_PAR  );	/*in, out, var*/

	void   shift_origin_inplace( double * , ENV_PAR  );				/*in-out*/
	void unshift_origin_inplace( double * , ENV_PAR  );				/*in-out*/

	void   shift_origin_complex_inplace( fftw_complex *, ENV_PAR  );		/*in-out*/
	void unshift_origin_complex_inplace( fftw_complex *, ENV_PAR  );		/*in-out*/

#endif

