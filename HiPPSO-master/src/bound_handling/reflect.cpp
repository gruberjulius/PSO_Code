/**
* @file   bound_handling/reflect.cpp
* @author Alexander Raß (alexander.rass@fau.de)
* @date   September, 2015
* @brief  This file contains the reflect bound handling strategy.
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

#include "bound_handling/reflect.h"

#include "function/function.h"
#include "general/configuration.h"
#include "arbitrary_precision_calculation/operations.h"
#include "general/particle.h"
#include "velocity_adjustment/velocity_adjustment.h"

namespace highprecisionpso {

void BoundHandlingReflect::SetParticleUpdate(Particle * p){
	std::vector<mpf_t*> oldPos = p->GetPosition();
	std::vector<mpf_t*> vel = p->GetVelocity();
	std::vector<mpf_t*> newPos = arbitraryprecisioncalculation::vectoroperations::Add(oldPos, vel);
	std::vector<mpf_t*> low_position = configuration::g_function->GetLowerSearchSpaceBound();
	std::vector<mpf_t*> high_position = configuration::g_function->GetUpperSearchSpaceBound();
	std::vector<bool> modifiedDimensions(newPos.size(), false);
	for(unsigned int d = 0; d < newPos.size(); d++){
		while(arbitraryprecisioncalculation::mpftoperations::Compare(newPos[d], low_position[d]) < 0 || arbitraryprecisioncalculation::mpftoperations::Compare(newPos[d], high_position[d]) > 0){

			modifiedDimensions[d] = true;

			mpf_t* bound = (mpf_t*) NULL;
			if(arbitraryprecisioncalculation::mpftoperations::Compare(newPos[d], low_position[d]) < 0) {
				bound = low_position[d];
			} else {
				bound = high_position[d];
			}
			// It is not allowed to free bound because it is not a cloned value.
			mpf_t* dif = arbitraryprecisioncalculation::mpftoperations::Subtract(newPos[d], bound);
			arbitraryprecisioncalculation::mpftoperations::ReleaseValue(newPos[d]);
			newPos[d] = arbitraryprecisioncalculation::mpftoperations::Subtract(bound, dif);
			arbitraryprecisioncalculation::mpftoperations::ReleaseValue(dif);
		}
	}
	p->SetPosition(newPos);
	configuration::g_velocity_adjustment->AdjustVelocity(p, modifiedDimensions, oldPos);
	arbitraryprecisioncalculation::vectoroperations::ReleaseValues(low_position);
	arbitraryprecisioncalculation::vectoroperations::ReleaseValues(high_position);
	arbitraryprecisioncalculation::vectoroperations::ReleaseValues(oldPos);
	arbitraryprecisioncalculation::vectoroperations::ReleaseValues(vel);
	arbitraryprecisioncalculation::vectoroperations::ReleaseValues(newPos);
}

std::string BoundHandlingReflect::GetName(){
	return "Ref";
}

} // namespace highprecisionpso
