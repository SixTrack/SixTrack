#ifndef FLUKAIO_FORTRANFLUKAIO_H__
#define FLUKAIO_FORTRANFLUKAIO_H__

#include <stdint.h>

#ifdef __cplusplus
extern "C" {
#endif

/**
 * Initialize FlukaIO **Mandatory**
 */
int ntinit_();

/**
 * wait and return number of particles to handle
 * @param cid Connection ID
 * @return number of particles to handle
 */
int ntnpart_(const int *cid);

/**
 * wait and return the Brho nominal
 * @param cid Connection ID
 * @return -1 on error 0 on ok
 */
int ntbrho_(const int *cid, double *brho);

/**
 * Receive message
 * All parameters are output.
 * @param cid Connection ID
 * @param type Message type as defined in Message.h
 * @param id Particle id
 * @param gen Generation
 * @param wgt Particle statistical weight
 * @param x
 * @param y
 * @param z
 * @param tx
 * @param ty
 * @param tz
 * @param aa
 * @param zz
 * @param m
 * @param pc
 * @return message size if ok, -1 if error (disconnect, timeout, etc.)
 */
int ntrecv_(
		const int *cid,
		uint8_t *type,
		uint32_t *id, uint32_t *gen,
		double *wgt,
		double  *x, double  *y, double *z,
		double *tx, double *ty, double *tz,
		uint16_t *aa, uint16_t *zz,
		double *m, double *pc,
		double *t);

/**
 * Wait message
 * All params are output. Blocks process until a message is read
 * @param cid Connection ID
 * @param type Message type as defined in Message.h
 * @param id Particle id
 * @param gen Generation
 * @param wgt Particle statistical weight
 * @param x
 * @param y
 * @param z
 * @param tx
 * @param ty
 * @param tz
 * @param aa
 * @param zz
 * @param m
 * @param pc
 * @return message size if ok, -1 if error (disconnect, timeout, etc.)
 */
int ntwait_(
		const int *cid,
		uint8_t *type,
		uint32_t *id, uint32_t *gen,
		double *wgt,
		double  *x, double  *y, double *z,
		double *tx, double *ty, double *tz,
		uint16_t *aa, uint16_t *zz,
		double *m, double *pc,
		double *t);

/**
 * Send particle
 * @param cid Connection ID
 * @param id Particle if
 * @param wgt Particle statistical weight
 * @param x
 * @param y
 * @param z
 * @param tx
 * @param ty
 * @param tz
 * @param aa
 * @param zz
 * @param m
 * @param pc
 * @return message size if ok, -1 if error
 */
int ntsendp_(
		const int *cid,
		const uint32_t *id, const uint32_t *gen,
		const double *wgt,
		const double  *x,   const double  *y, const double *z,
		const double *tx,   const double *ty, const double *tz,
		const uint16_t *aa, const uint16_t *zz,
		const double *m, const double *pc,
		const double *t);

/**
 * Send End of Turn
 * @param cid Connection ID
 * @return message size if ok, -1 if error
 */
int ntsendeob_(const int *cid);

/**
 * Send Insertion Point
 * @param cid Connection ID
 * @param turn Turn number of the current turn
 * @param ipt Insertion Point number
 * @return message size if ok, -1 if error
 */
int ntsendipt_(const int *cid, const uint32_t *turn, const uint16_t *ipt);

/**
 * Send End of Computation
 * @param cid Connection ID
 * @return message size if ok, -1 if error
 */
int ntsendeoc_(const int *cid);

/**
 * Send maximum number of particles
 * @param npart number of particles
 * @return message size if ok, -1 if error
 */
int ntsendnpart_(const int *cid, const uint32_t *npart);

/**
 * Send nominal Brho
 * @param brhono nominal Brho
 * @return message size if ok, -1 if error
 */
int ntsendnbrhono_(const int *cid, const double *brhono);

int ntbrho_(const int *cid, double *brhono);

/**
 * Set read timeout
 * @param cid Connection ID
 * @param seconds
 * @return -1 if error
 */
int ntrtimeout_(const int *cid, const int *seconds);

/**
 * Set write timeout
 * @param cid Connection ID
 * @param seconds
 * @return -1 if error
 */
int ntwtimeout_(const int *cid, const int *seconds);

/* ***************************************
 * Server management
 * ***************************************/

/**
 * Create a server instance and retunr its identifier
 * @return server id or -1 ir error
 */
int ntserver_();

/**
 * Start server and listen in port
 * @param port port where to listen, (randomly assigned by system if 0)
 * @return assigned port if ok, -1 if error
 */
int ntstart_(const int * serverid, int *port);

/**
 * Blocks until incoming connection (as server)
 * @return connection id, -1 if error
 */
int ntaccept_(const int *serverid);

/**
 * Shutdown server
 * Releases server resources
 * @return 0 if ok, -1 if error
 */
int ntshdwn_(const int *serverid);

/**
 * Get server Port
 * @return Port number where the server is listening
 */
int ntgetport_(const int *serverid);

/**
 * Connect to server
 * @param host
 * @param port
 * @return connection id if success, -1 if error
 */
int ntconnect_(char *host, int *port, const long hostlen);

/**
 * Finalize connection
 * @param cid connection id
 * @return 0 if success, -1 if error
 */
int ntend_(const int *cid);

char * create_fortran_string(const char *str, const long str_len);

#ifdef __cplusplus
}
#endif

#endif
