#ifndef FLUKAIO_PARTICLE_INFO_H__
#define FLUKAIO_PARTICLE_INFO_H__

#include <stdint.h>

#ifdef __cplusplus
extern "C" {
#endif

/**
 * Stores particle data as sent through the network 
 * Always using Fluka units
 *
 * Fields:
 *
 * - Particle identification:
 *    - id: Particle id
 *    - gen: Particle generation
 * - Particle info:
 *    - x: Position in x (cm.)
 *    - y: Position in y (cm.)
 *    - z: Position in z (cm.)
 *    - tx: Director cosine in x
 *    - ty: Director cosine in y
 *    - tz: Director cosine in z
 *    - m: Rest mass (GeV/c^2)
 *    - p: Momentum (GeV/c)
 *    - aa: Mass number (total number of protons + neutrons)
 *    - zz: Ion charge (number of protons)
 *    - t: Time
 * - Particle metadata:
 *    - Weight: statistical weight (goes from 0.0 to 1.0)
 *
 */

#pragma pack(1)
typedef struct {
	/* Particle identification */
	uint32_t id;
	uint32_t gen;

	/* Particle attributtes */
	uint16_t aa;
	uint16_t zz;
	double x;
	double y;
	double z;
	double tx;
	double ty;
	double tz;
	double m;
	double pc;
	double t;

	/* Particle metadata */
	double weight;
} particle_info_t;
#pragma pack()

#ifdef __cplusplus
}
#endif

#endif
