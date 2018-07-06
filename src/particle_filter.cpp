/*
 * particle_filter.cpp
 *
 *  Created on: Dec 12, 2016
 *      Author: Tiffany Huang
 */

#include <random>
#include <algorithm>
#include <iostream>
#include <numeric>
#include <math.h> 
#include <iostream>
#include <sstream>
#include <string>
#include <iterator>

#include "particle_filter.h"

#define EPS 0.0001

using namespace std;

void ParticleFilter::init(double x, double y, double theta, double std[]) {
	// TODO: Set the number of particles. Initialize all particles to first position (based on estimates of 
	//   x, y, theta and their uncertainties from GPS) and all weights to 1. 
	// Add random Gaussian noise to each particle.
	// NOTE: Consult particle_filter.h for more information about this method (and others in this file).

	num_particles = 100;
	default_random_engine gen;
	
	// set distribution for input data. 
	// 1. x 
	// 2. y
	// 3. theta
	normal_distribution<double> dist_x(x, std[0]);
	normal_distribution<double> dist_y(y, std[1]);
	normal_distribution<double> dist_theta(theta, std[2]);


	for (int i = 0; i < num_particles; i++) {

		Particle curr_particle;
		curr_particle.id = i;

		curr_particle.x = dist_x(gen);
		curr_particle.y = dist_y(gen);
		curr_particle.theta = dist_theta(gen);
		curr_particle.weight = 1.0;

		particles.push_back(curr_particle);
		weights.push_back(curr_particle.weight);
	}
	is_initialized = true;

}

void ParticleFilter::prediction(double delta_t, double std_pos[], double velocity, double yaw_rate) {
	// TODO: Add measurements to each particle and add random Gaussian noise.
	// NOTE: When adding noise you may find std::normal_distribution and std::default_random_engine useful.
	//  http://en.cppreference.com/w/cpp/numeric/random/normal_distribution
	//  http://www.cplusplus.com/reference/random/default_random_engine/

	default_random_engine gen;

	normal_distribution<double> dist_x(0, std_pos[0]);
	normal_distribution<double> dist_y(0, std_pos[1]);
	normal_distribution<double> dist_theta(0, std_pos[2]);


	int i;
	for (i = 0; i < num_particles; i++) {

		double particle_theta = particles[i].theta;

		//Instead of a hard check of 0, adding a check for very low value of yaw_rate
		if (fabs(yaw_rate) < EPS) {

			particles[i].x += velocity * cos(particle_theta) * delta_t;
			particles[i].y += velocity * sin(particle_theta) * delta_t;

		} else {

			particles[i].x += (velocity/yaw_rate) * (sin(particle_theta + (yaw_rate * delta_t)) - sin(particle_theta));
			particles[i].y += (velocity/yaw_rate) * (cos(particle_theta) - cos(particle_theta + (yaw_rate * delta_t)));
			particles[i].theta += (yaw_rate * delta_t);

		}

		// add noise based on entered standard deviation.
		particles[i].x += dist_x(gen);
		particles[i].y += dist_y(gen);
		particles[i].theta += dist_theta(gen);

	}



}

void ParticleFilter::dataAssociation(std::vector<LandmarkObs> predicted, std::vector<LandmarkObs>& observations) {
	// TODO: Find the predicted measurement that is closest to each observed measurement and assign the 
	//   observed measurement to this particular landmark.
	// NOTE: this method will NOT be called by the grading code. But you will probably find it useful to 
	//   implement this method and use it as a helper during the updateWeights phase.

	
	int predicted_size = predicted.size();
	int observations_size = observations.size();

	for (int i=0;i<observations_size;i++ ) {

		// init minimum distance to maximum possible
		double min_dist = numeric_limits<double>::max();

		// pickup observations index		
		double o_x = observations[i].x;
		double o_y = observations[i].y;
		
		int mapid = -1;

		for (int j=0;j<predicted_size;j++) {

			double pred_x = predicted[j].x;
			double pred_y = predicted[j].y;
			 
			double distance = dist( o_x,o_y,pred_x,pred_y );

			if (distance < min_dist) {
				min_dist = distance;
				mapid = predicted[j].id; // find most nearest distance mapid
			}
		}

		observations[i].id = mapid;
	}
	

}



void ParticleFilter::updateWeights(double sensor_range, double std_landmark[], 
		const std::vector<LandmarkObs> &observations, const Map &map_landmarks) {
	// TODO: Update the weights of each particle using a mult-variate Gaussian distribution. You can read
	//   more about this distribution here: https://en.wikipedia.org/wiki/Multivariate_normal_distribution
	// NOTE: The observations are given in the VEHICLE'S coordinate system. Your particles are located
	//   according to the MAP'S coordinate system. You will need to transform between the two systems.
	//   Keep in mind that this transformation requires both rotation AND translation (but no scaling).
	//   The following is a good resource for the theory:
	//   https://www.willamette.edu/~gorr/classes/GeneralGraphics/Transforms/transforms2d.htm
	//   and the following is a good resource for the actual equation to implement (look at equation 
	//   3.33
	//   http://planning.cs.uiuc.edu/node99.html


	// 
	//double weight_normalizer = 0.0;
	// loop number of particles
	for (int i=0;i<num_particles;i++) {

		//
		// save main particles member values
		//
		double p_x = particles[i].x;
		double p_y = particles[i].y;
		double p_theta = particles[i].theta;


		// step 1 
		// observations are transformed from vehicles co-ordinaties to map co-ordinates 		
		std::vector<LandmarkObs> transformed_observations;
		for (int j=0;j<observations.size();j++) {
			LandmarkObs trans_obs;
			trans_obs.x = (cos(p_theta) * observations[j].x - sin(p_theta) * observations[j].y) + p_x;
			trans_obs.y = (sin(p_theta) * observations[j].x + cos(p_theta) * observations[j].y) + p_y;
			trans_obs.id = observations[j].id;
			// saved to new transformed matrix 
			transformed_observations.push_back( trans_obs  );
		}


		//
		// step 2
		// find target landmark dx ^ 2 + dy ^ 2 less than sensor_range ^ 2
		// set most nearest landmark point to predicted_landmarks
		// 
		std::vector<LandmarkObs> predicted_Landmarks;
		int landmark_size = map_landmarks.landmark_list.size();		
		for (int j=0;j<landmark_size;j++) {

			double lm_x = map_landmarks.landmark_list[j].x_f;
			double lm_y = map_landmarks.landmark_list[j].y_f;
			int lm_id = map_landmarks.landmark_list[j].id_i;

			//double dx = p_x - lm_x;
			//double dy = p_y - lm_y;

      		if ((fabs(( p_x - lm_x)) <= sensor_range) && (fabs( ( p_y - lm_y) ) <=  sensor_range) ) {
			//if ( (pow(dx,2) + pow(dy,2)) < pow(sensor_range,2 ) ) {
				// if distance is ranged inside sensor_range 
				// save landmark id x y into new vector..
				predicted_Landmarks.push_back( LandmarkObs{ lm_id,lm_x,lm_y }     );

			}
		}

		//

		//
		// step 3 
		// call data assosciations.
		dataAssociation(  predicted_Landmarks, transformed_observations  );

		// step 4 
		// calculate weights
		// 1. targetLandmaks (within range)
		// 2. transformed observations (based on particles x y theta)
		particles[i].weight = 1.0;

		double x_sigma = std_landmark[0];
		double y_sigma = std_landmark[1];
		double x_sigma2 = pow(x_sigma,2);
		double y_sigma2 = pow(y_sigma,2);
		double multiv1 = (1.0 / ( 2 * M_PI * x_sigma * y_sigma  ));		
		
		for (int k=0;k<transformed_observations.size();k++) {

			double trans_o_x = transformed_observations[k].x;
			double trans_o_y = transformed_observations[k].y;
			int trans_o_id = transformed_observations[k].id;

			double probs = 1.0;
			for (int l=0; l < predicted_Landmarks.size(); l++) {

				double pred_x = predicted_Landmarks[l].x;
				double pred_y = predicted_Landmarks[l].y;
				int pred_id = predicted_Landmarks[l].id;

				if (trans_o_id == pred_id) {

					double dx2 = pow( (trans_o_x - pred_x), 2 );
					double dxCalc = ( dx2 /  (2. * x_sigma2 ) );

					double dy2 = pow( (trans_o_y - pred_y), 2 );
					double dyCalc = ( dy2 / ( 2. * y_sigma2 ) );

					probs = multiv1 * exp( -1.0 * ( dxCalc + dyCalc ) );
					particles[i].weight *= probs;

					//std::cout << "calculate  (idx):" << k << ":" << l << "probs" << probs << std::endl;
				}

			}
		}
    	//weight_normalizer += particles[i].weight;

	// end of particles loop..
	}



}

void ParticleFilter::resample() {
	// TODO: Resample particles with replacement with probability proportional to their weight. 
	// NOTE: You may find std::discrete_distribution helpful here.
	//   http://en.cppreference.com/w/cpp/numeric/random/discrete_distribution



	/* copy weights from particle weights */
	for (int i = 0; i < particles.size(); i++) {
		//particles[i].weight /= weight_normalizer;
		weights[i] = particles[i].weight;
	}

	double max_weight = *std::max_element(weights.begin(), weights.end());
	double max_weight2 = 2.0 * max_weight;
	double beta = 0.0f;

	// build uniform random 
	// 1 weight random
	// 2 index 
	uniform_real_distribution<double> random_weight(0.0, max_weight2);
	uniform_int_distribution<int> random_idx(0,num_particles-1);

	// Create a generator to be used for generating random particle index and beta value
	default_random_engine gen;

	std::vector<Particle> resampledParticles;
	int index = random_idx(gen);

	// wheel process
	for (int i=0; i < num_particles; i++) {

		beta += random_weight(gen);
		while (beta > weights[index]) {
			
			beta -= weights[index];
			index = (index + 1) % num_particles;

		}
		resampledParticles.push_back( particles[index] );

	}
	particles = resampledParticles;
	
}

Particle ParticleFilter::SetAssociations(Particle& particle, const std::vector<int>& associations, 
                                     const std::vector<double>& sense_x, const std::vector<double>& sense_y)
{
    //particle: the particle to assign each listed association, and association's (x,y) world coordinates mapping to
    // associations: The landmark id that goes along with each listed association
    // sense_x: the associations x mapping already converted to world coordinates
    // sense_y: the associations y mapping already converted to world coordinates

    particle.associations= associations;
    particle.sense_x = sense_x;
    particle.sense_y = sense_y;

	return particle;
}

string ParticleFilter::getAssociations(Particle best)
{
	vector<int> v = best.associations;
	stringstream ss;
    copy( v.begin(), v.end(), ostream_iterator<int>(ss, " "));
    string s = ss.str();
    s = s.substr(0, s.length()-1);  // get rid of the trailing space
    return s;
}
string ParticleFilter::getSenseX(Particle best)
{
	vector<double> v = best.sense_x;
	stringstream ss;
    copy( v.begin(), v.end(), ostream_iterator<float>(ss, " "));
    string s = ss.str();
    s = s.substr(0, s.length()-1);  // get rid of the trailing space
    return s;
}
string ParticleFilter::getSenseY(Particle best)
{
	vector<double> v = best.sense_y;
	stringstream ss;
    copy( v.begin(), v.end(), ostream_iterator<float>(ss, " "));
    string s = ss.str();
    s = s.substr(0, s.length()-1);  // get rid of the trailing space
    return s;
}
