/**
* @file   bound_handling/torus.h
* @author Alexander Raß (alexander.rass@fau.de)
* @date   September, 2015
* @brief  This file contains the torus bound handling strategy.
*
* @copyright
* This project is released under the MIT License (MIT).
*
* @copyright
* The MIT License (MIT)
*
* @copyright
* Copyright (c) 2016 by Friedrich-Alexander-Universität Erlangen-Nürnberg and
* Alexander Raß
*
* @copyright
* Permission is hereby granted, free of charge, to any person obtaining a copy
* of this software and associated documentation files (the "Software"), to deal
* in the Software without restriction, including without limitation the rights
* to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
* copies of the Software, and to permit persons to whom the Software is
* furnished to do so, subject to the following conditions:
*
* @copyright
* The above copyright notice and this permission notice shall be included in
* all copies or substantial portions of the Software.
*
* @copyright
* THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
* IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
* FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
* AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
* LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
* OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
* SOFTWARE.
*
*/

#ifndef HIGH_PRECISION_PSO_BOUND_HANDLING_TORUS_H_
#define HIGH_PRECISION_PSO_BOUND_HANDLING_TORUS_H_

#include "bound_handling/bound_handling.h"

namespace highprecisionpso {

/**
* @brief This class implements the bound handling strategy torus for the particle swarm optimization.
*
* If a particle would leave the search space on one side, it will enter the search space on the other side.
* The movement is similar to modular arithmetic.
*/
class BoundHandlingTorus : public BoundHandling{
public:
	/**
	* @brief Updates the position of the particle p according to the current position and the new (already calculated) velocity.
	*
	* The new position is calculated by adding the velocity to the position.
	* If the new position would be outside of the search space bounds then the bound handling strategy torus is activated.
	* If a particle would leave the search space on one side, it will enter the search space on the other side.
	* The movement is similar to modular arithmetic.
	* No velocity adjustment is executed here.
	* 
	* @param p The particle, which should be updated.
	*/
	void SetParticleUpdate(Particle * p);
	std::string GetName();
	/**
	* @brief Calculates the difference vector between the position and the aim.
	*
	* As the search space repeats, the aim can be reached in positive direction as well as in negative direction.
	* The value with smaller absolute value counts as the difference in that dimension.
	*
	* @param position The vector, which describes the position.
	* @param aim The vector, which describes the aim.
	*
	* @return The difference from the position to the aim.
	*/
	std::vector<mpf_t*> GetDirectionVector(const std::vector<mpf_t*> & position, const std::vector<mpf_t*> & aim);
};

} // namespace highprecisionpso

#endif /* HIGH_PRECISION_PSO_BOUND_HANDLING_TORUS_H_ */
