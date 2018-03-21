#ifndef FLUKAIO_MESSAGE_H__
#define FLUKAIO_MESSAGE_H__

#include <stdint.h>
#include <stdlib.h>
#include "ParticleInfo.h"

#ifdef __cplusplus
extern "C" {
#endif

/** Message type type */
typedef uint8_t flukaio_message_type_t;

#pragma pack(1)
typedef struct { // Insertion point message
	uint32_t turn;
	uint16_t ipt;
} flukaio_ipt_data_t;

typedef struct { // Handshake
	uint16_t major; // version
	uint16_t minor;
	uint32_t key; // random key
} flukaio_hsk_data_t;

typedef struct { // integer variable passing
	uint16_t id;	// variable id
	int32_t	 value;
} flukaio_int_data_t;

typedef struct { // double variable passing
	uint16_t id;	// variable id
	double	 value;
} flukaio_dble_data_t;
#pragma pack()

enum flukaio_int_id_e {
	N_NPART  = 0x0
};

enum flukaio_dble_id_e {
	N_BRHONO = 0x0
};

// Message lengths
#define MSG_HEADER_LEN  (sizeof(uint16_t)+sizeof(flukaio_message_type_t))

/** Message types */
enum flukaio_message_type_e {
	N_ERR   = 0x00,
	N_PART  = 0x01,
	N_EOB   = 0x02,
	N_END   = 0x03,
	N_CONF  = 0x04,
	N_IPT   = 0x05,
	N_HSK   = 0x06,
	N_INT   = 0x07,
	N_DBLE  = 0x08,
	N_OTHER = 0xff
};

/**
 * Network message
 *
 * With uint16_t size the maximum size for data is 2^16-3 (3 is header length)
 */
#pragma pack(1)
typedef struct {
	uint16_t size; /**< Total size in bytes of the message (including headers) */
	flukaio_message_type_t type;  /**< Type of message */

	// data field, depending on the size, can contain a particle, an ipt or hsk
	union _data {
		particle_info_t       particle;
		flukaio_ipt_data_t    ipt;
		flukaio_hsk_data_t    hsk;
		flukaio_int_data_t    varint;
		flukaio_dble_data_t   vardouble;
	} data;

} flukaio_message_t;
#pragma pack()

#ifdef __cplusplus
}
#endif

#endif
