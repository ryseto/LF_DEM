//
//  Averager.h
//  LF_DEM
//
//  Created by Romain Mari on 02/04/15.
//  Copyright (c) 2015 Ryohei Seto and Romain Mari. All rights reserved.
//

/**
 \class Averager
 \brief Small class computing and storing observable averages with an exponential memory kernel.

 The average x_avg of quantity x is defined as:
 \f$x_{avg}(\gamma_n) = C_n \sum_{i=1}^n x(\gamma_i)e^{(\gamma_n-\gamma_i)/\gamma_{avg}} \f$

 \author Ryohei Seto
 \author Romain Mari
 */


#ifndef __LF_DEM__Averager__
#define __LF_DEM__Averager__

template <typename XType>
class Averager{
 private:
	XType x_avg;
	double previous_kernel_norm;
	double relaxation_time;
 public:
	void update(XType x, double deltagamma){
		/**	  
		 * \brief Updates the average value x_avg with an exponential memory kernel
		 * 
		 * The average x_avg is defined as:
		 * \f$x_{avg}(\gamma_n) = C_n \sum_{i=1}^n x(\gamma_i)e^{(\gamma_n-\gamma_i)/\gamma_{avg}} \f$
		 *
		 * \param  x : \f$x_n \f$ 
		 * \param  x_avg : \f$x_{avg}(\gamma_{n-1}) \f$ (at strain \f$n-1\f$)
		 * \param  deltagamma : \f$\gamma_n-\gamma_{n-1} \f$.
		 * 
		 * On return, x_avg contains the value of \f$x_{avg}(\gamma_{n}) \f$ (at strain \f$n\f$)
		 *
		 * This method uses the following recursions:
		 * - \f$ C_n = \frac{C_{n-1}}{(\gamma_n-\gamma_{n-1})C_{n-1}} + e^{(\gamma_n-\gamma_{n-1})/\gamma_{avg}} \f$
		 * - \f$ x_{avg}(\gamma_n) = C_{n}\left[ (\gamma_n-\gamma_{n-1}) + e^{(\gamma_n-\gamma_{n-1})/\gamma_{avg}} x_{avg}(\gamma_{n-1})/C_{n-1} \right] \f$
		 */
		
		double etn = exp(-deltagamma/relaxation_time);
		double inv_kernel_norm = previous_kernel_norm/(deltagamma*previous_kernel_norm+etn);
		if (previous_kernel_norm == 1) { // = it is the first iteration
			inv_kernel_norm=1/deltagamma;
		}
		
		x_avg = inv_kernel_norm*( x*deltagamma + etn*x_avg/previous_kernel_norm );
		
		previous_kernel_norm = inv_kernel_norm;
	};
	
 
 Averager(double rtime) :
	x_avg(0),
		previous_kernel_norm(1),
		relaxation_time(rtime)
	{};
	~Averager();
	
	XType get(){ return x_avg;};
};
#endif /* defined(__LF_DEM__Averager__) */
