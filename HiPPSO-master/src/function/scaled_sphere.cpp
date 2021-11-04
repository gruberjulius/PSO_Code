/**
* @file   function/scaled_sphere.cpp
* @author Alexander Raß (alexander.rass@fau.de)
* @date   January, 2014
* @brief  This file contains the description of the scaled sphere functions.
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

#include <arbitrary_precision_calculation/random_number_generator.h>
#include "function/scaled_sphere.h"

#include <cmath>
#include <sstream>

#include "general/check_condition.h"
#include "arbitrary_precision_calculation/operations.h"
#include "arbitrary_precision_calculation/parse.h"

namespace highprecisionpso {

// ScaledSphere ScaledSphere ScaledSphere ScaledSphere ScaledSphere ScaledSphere ScaledSphere ScaledSphere
ScaledSphere::ScaledSphere(double maxScale) : max_scale_ (maxScale) {
	scales_cached_.clear();
	scales_cached_precision_ = 1;
}

ScaledSphere::~ScaledSphere(){
	arbitraryprecisioncalculation::mpftoperations::ChangeNumberOfMpftValuesCached( -((int)(scales_cached_.size())) );
	arbitraryprecisioncalculation::vectoroperations::ReleaseValues(scales_cached_);
	scales_cached_.clear();
}

mpf_t* ScaledSphere::Eval(const std::vector<mpf_t*> & vec) {
	unsigned int D = vec.size();
	if(scales_cached_.size() != D || scales_cached_precision_ != mpf_get_default_prec()) {
		scales_cached_precision_ = mpf_get_default_prec();
		arbitraryprecisioncalculation::mpftoperations::ChangeNumberOfMpftValuesCached( -((int)(scales_cached_.size())) );
		arbitraryprecisioncalculation::vectoroperations::ReleaseValues(scales_cached_);
		scales_cached_.clear();
		scales_cached_.push_back(arbitraryprecisioncalculation::mpftoperations::ToMpft(1.0));
		mpf_t* divisor = arbitraryprecisioncalculation::mpftoperations::ToMpft(std::max(1, ((int) D) - 1));
		mpf_t* one = arbitraryprecisioncalculation::mpftoperations::ToMpft(1);
		mpf_t* exponent = arbitraryprecisioncalculation::mpftoperations::Divide(one, divisor);
		mpf_t* base = arbitraryprecisioncalculation::mpftoperations::ToMpft(max_scale_);
		mpf_t* multiplier = arbitraryprecisioncalculation::mpftoperations::Pow(base, exponent);
		while(scales_cached_.size() < D){
			scales_cached_.push_back(arbitraryprecisioncalculation::mpftoperations::Multiply(scales_cached_[scales_cached_.size() - 1], multiplier));
		}
		arbitraryprecisioncalculation::mpftoperations::ReleaseValue(divisor);
		arbitraryprecisioncalculation::mpftoperations::ReleaseValue(one);
		arbitraryprecisioncalculation::mpftoperations::ReleaseValue(exponent);
		arbitraryprecisioncalculation::mpftoperations::ReleaseValue(base);
		arbitraryprecisioncalculation::mpftoperations::ReleaseValue(multiplier);
		arbitraryprecisioncalculation::mpftoperations::ChangeNumberOfMpftValuesCached( scales_cached_.size() );
	}
	std::vector<mpf_t*> vec2 = arbitraryprecisioncalculation::vectoroperations::Multiply(vec, vec);
	std::vector<mpf_t*> sumUpValues = arbitraryprecisioncalculation::vectoroperations::Multiply(vec2, scales_cached_);
	mpf_t* res = arbitraryprecisioncalculation::vectoroperations::Add(sumUpValues);
	arbitraryprecisioncalculation::vectoroperations::ReleaseValues(vec2);
	arbitraryprecisioncalculation::vectoroperations::ReleaseValues(sumUpValues);
	return res;
}

std::string ScaledSphere::GetName(){
	std::ostringstream os;
	os << "Sc2NV" << max_scale_;
	return os.str();
}

mpf_t* ScaledSphere::DistanceTo1DLocalOptimum(const std::vector<mpf_t*> & pos, int d){
	return arbitraryprecisioncalculation::mpftoperations::Abs(pos[d]);
}
// ScaledSphereFixed ScaledSphereFixed ScaledSphereFixed ScaledSphereFixed ScaledSphereFixed ScaledSphereFixed
ScaledSphereFixed::ScaledSphereFixed():ScaledSphere(1000000.0){}

// ScaledHadamardRotatedSphere ScaledHadamardRotatedSphere ScaledHadamardRotatedSphere ScaledHadamardRotatedSphere
ScaledHadamardRotatedSphere::ScaledHadamardRotatedSphere(double max_scale) : max_scale_ (max_scale){
	matrix_cached_.clear();
	matrix_cached_precision_ = 1;
}

ScaledHadamardRotatedSphere::~ScaledHadamardRotatedSphere(){
	for(unsigned int i = 0; i < matrix_cached_.size(); i++){
		arbitraryprecisioncalculation::mpftoperations::ChangeNumberOfMpftValuesCached( - ((int)(matrix_cached_[i].size())) );
		arbitraryprecisioncalculation::vectoroperations::ReleaseValues(matrix_cached_[i]);
		matrix_cached_[i].clear();
	}
}

void ScaledHadamardRotatedSphere::InitMatrixCached(const std::vector<mpf_t*> & vec) {
	unsigned int D = vec.size();
	if(matrix_cached_.size() != D || matrix_cached_precision_ != mpf_get_default_prec()){
		matrix_cached_precision_ = mpf_get_default_prec();
		for(unsigned int i = 0; i < matrix_cached_.size(); i++){
			arbitraryprecisioncalculation::mpftoperations::ChangeNumberOfMpftValuesCached( -((int)(matrix_cached_[i].size())) );
			arbitraryprecisioncalculation::vectoroperations::ReleaseValues(matrix_cached_[i]);
			matrix_cached_[i].clear();
		}
		mpf_t* divisor = arbitraryprecisioncalculation::mpftoperations::ToMpft(2 * std::max(1, ((int) D) - 1));
		mpf_t* one = arbitraryprecisioncalculation::mpftoperations::ToMpft(1);
		mpf_t* exponent = arbitraryprecisioncalculation::mpftoperations::Divide(one, divisor);
		mpf_t* base = arbitraryprecisioncalculation::mpftoperations::ToMpft(max_scale_);
		mpf_t* multiplier = arbitraryprecisioncalculation::mpftoperations::Pow(base, exponent);
		arbitraryprecisioncalculation::mpftoperations::ReleaseValue(divisor);
		arbitraryprecisioncalculation::mpftoperations::ReleaseValue(one);
		arbitraryprecisioncalculation::mpftoperations::ReleaseValue(exponent);
		arbitraryprecisioncalculation::mpftoperations::ReleaseValue(base);
		std::vector<std::vector<mpf_t*> > tmp(D, std::vector<mpf_t*>(D, (mpf_t*) NULL));
		for(unsigned int i = 0; i < D; i++){
			for(unsigned int j = 0; j < D; j++){
				double neg = 1.0;
				int tmp_and = i&j;
				while(tmp_and){
					neg = -neg;
					tmp_and -= tmp_and&-tmp_and;
				}
				tmp[j][i] = arbitraryprecisioncalculation::mpftoperations::Pow(multiplier, (int)i);
				if(neg < 0){
					mpf_t* temp = arbitraryprecisioncalculation::mpftoperations::Negate(tmp[j][i]);
					std::swap(temp, tmp[j][i]);
					arbitraryprecisioncalculation::mpftoperations::ReleaseValue(temp);
				}
			}
		}
		arbitraryprecisioncalculation::mpftoperations::ReleaseValue(multiplier);
		matrix_cached_ = std::vector<std::vector<mpf_t*> >(D, std::vector<mpf_t*>(D, (mpf_t*) NULL));
		for(unsigned int i = 0; i < D; i++){
			for(unsigned int j = 0; j < D; j++){
				matrix_cached_[i][j] = arbitraryprecisioncalculation::mpftoperations::ToMpft(0.0);
				for(unsigned int k = 0; k < D; k++){
					mpf_t* next1 = arbitraryprecisioncalculation::mpftoperations::Multiply(tmp[i][k], tmp[j][k]);
					mpf_t* next2 = arbitraryprecisioncalculation::mpftoperations::Add(matrix_cached_[i][j], next1);
					arbitraryprecisioncalculation::mpftoperations::ReleaseValue(matrix_cached_[i][j]);
					arbitraryprecisioncalculation::mpftoperations::ReleaseValue(next1);
					matrix_cached_[i][j] = next2;
				}
			}
		}
		arbitraryprecisioncalculation::mpftoperations::ChangeNumberOfMpftValuesCached ( D * D );
		for(unsigned int i = 0; i < D; i++){
			arbitraryprecisioncalculation::vectoroperations::ReleaseValues(tmp[i]);
		}
	}
}

mpf_t* ScaledHadamardRotatedSphere::Eval(const std::vector<mpf_t*> & vec) {
	InitMatrixCached(vec);
	unsigned int D = vec.size();
	mpf_t* res = arbitraryprecisioncalculation::mpftoperations::ToMpft(0.0);
	for(unsigned int i = 0; i < D; i++){
		for(unsigned int j = 0;j < D; j++){
			mpf_t* cur = arbitraryprecisioncalculation::mpftoperations::Multiply(matrix_cached_[i][j], vec[i]);
			mpf_t* next = arbitraryprecisioncalculation::mpftoperations::Multiply(cur, vec[j]);
			arbitraryprecisioncalculation::mpftoperations::ReleaseValue(cur);
			cur = next;
			next = arbitraryprecisioncalculation::mpftoperations::Add(cur, res);
			arbitraryprecisioncalculation::mpftoperations::ReleaseValue(cur);
			arbitraryprecisioncalculation::mpftoperations::ReleaseValue(res);
			res = next;
		}
	}
	return res;
}

std::string ScaledHadamardRotatedSphere::GetName(){
	std::ostringstream os;
	os << "ScHadRot2NV" << max_scale_;
	return os.str();
}

mpf_t* ScaledHadamardRotatedSphere::DistanceTo1DLocalOptimum(const std::vector<mpf_t*> & pos, int d){
	InitMatrixCached(pos);
	int D = pos.size();
	AssertCondition(D == (int)matrix_cached_.size() && matrix_cached_precision_ == mpf_get_default_prec(), "scaledandhadamardrotatedshpere: Initialization failed.");
	mpf_t* B = arbitraryprecisioncalculation::mpftoperations::ToMpft(0.0);
	for(int d2 = 0; d2 < D; d2++)if(d2 != d){
		mpf_t* f = matrix_cached_[d][d2];
		mpf_t* b = arbitraryprecisioncalculation::mpftoperations::Multiply(f,pos[d2]);
		mpf_t* newB = arbitraryprecisioncalculation::mpftoperations::Add(b,B);
		arbitraryprecisioncalculation::mpftoperations::ReleaseValue(b);
		arbitraryprecisioncalculation::mpftoperations::ReleaseValue(B);
		B = newB;
	}
	mpf_t* best = arbitraryprecisioncalculation::mpftoperations::Divide(B,matrix_cached_[d][d]);
	arbitraryprecisioncalculation::mpftoperations::ReleaseValue(B);
	mpf_t* res = arbitraryprecisioncalculation::mpftoperations::Add(best, pos[d]);
	arbitraryprecisioncalculation::mpftoperations::ReleaseValue(best);
	mpf_t* res_abs = arbitraryprecisioncalculation::mpftoperations::Abs(res);
	arbitraryprecisioncalculation::mpftoperations::ReleaseValue(res);
	return res_abs;
}

// ScaledHadamardRotatedSphereFixed ScaledHadamardRotatedSphereFixed ScaledHadamardRotatedSphereFixed
ScaledHadamardRotatedSphereFixed::ScaledHadamardRotatedSphereFixed():ScaledHadamardRotatedSphere(1000000.0){}


// ScaledSphereRand ScaledSphereRand ScaledSphereRand ScaledSphereRand ScaledSphereRand ScaledSphereRand ScaledSphereRand
ScaledSphereRand::ScaledSphereRand(std::vector<std::string> rng_description):rng_description_(rng_description){
	scales_cached_.clear();
	scales_cached_precision_ = 0;
}

ScaledSphereRand::~ScaledSphereRand(){
	arbitraryprecisioncalculation::mpftoperations::ChangeNumberOfMpftValuesCached( -((int)(scales_cached_.size())) );
	arbitraryprecisioncalculation::vectoroperations::ReleaseValues(scales_cached_);
	scales_cached_.clear();
}

mpf_t* ScaledSphereRand::Eval(const std::vector<mpf_t*> & vec) {
	unsigned int D = vec.size();
	if(scales_cached_.size() != D || scales_cached_precision_ != mpf_get_default_prec()){
		scales_cached_precision_ = mpf_get_default_prec();
		arbitraryprecisioncalculation::mpftoperations::ChangeNumberOfMpftValuesCached( -((int)(scales_cached_.size())) );
		arbitraryprecisioncalculation::vectoroperations::ReleaseValues(scales_cached_);
		scales_cached_.clear();
		arbitraryprecisioncalculation::RandomNumberGenerator* rng = NULL;
		{
			unsigned int parsed = 0;
			rng = arbitraryprecisioncalculation::parse::ParseRandomNumberGenerator(rng_description_, parsed);
		}
		while(scales_cached_.size() < D){
			scales_cached_.push_back(arbitraryprecisioncalculation::mpftoperations::ToMpft(std::max(1LL, rng->RandomLongLong())));
		}
		arbitraryprecisioncalculation::mpftoperations::ChangeNumberOfMpftValuesCached( scales_cached_.size() );
	}
	mpf_t* res = arbitraryprecisioncalculation::mpftoperations::ToMpft(0.0);
	for(unsigned int i = 0; i < D; i++){
		mpf_t* cur = arbitraryprecisioncalculation::mpftoperations::Multiply(scales_cached_[i], vec[i]);
		mpf_t* next = arbitraryprecisioncalculation::mpftoperations::Multiply(cur, vec[i]);
		arbitraryprecisioncalculation::mpftoperations::ReleaseValue(cur);
		cur = next;
		next = arbitraryprecisioncalculation::mpftoperations::Add(cur, res);
		arbitraryprecisioncalculation::mpftoperations::ReleaseValue(cur);
		arbitraryprecisioncalculation::mpftoperations::ReleaseValue(res);
		res = next;
	}
	return res;
}

std::string ScaledSphereRand::GetName(){
	std::ostringstream os;
	os << "Sc2NRand";
	return os.str();
}

mpf_t* ScaledSphereRand::DistanceTo1DLocalOptimum(const std::vector<mpf_t*> & pos, int d){
	return arbitraryprecisioncalculation::mpftoperations::Abs(pos[d]);
}

} // namespace highprecisionpso
