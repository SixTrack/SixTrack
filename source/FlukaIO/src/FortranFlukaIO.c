#include "FortranFlukaIO.h"

#include "FlukaIO.h"
#include "FlukaIOServer.h"
#include <string.h>

#include <stdio.h>

#define FLUKAIO_MAX_CONNECTIONS 15
#define FLUKAIO_MAX_SERVERS 5

/*
 * Fortran can't handle pointers therefore an array storing the pointers is created
 * and array indexes are given to fortran
 */
static flukaio_connection_t *connections[FLUKAIO_MAX_CONNECTIONS];
static flukaio_server_t *servers[FLUKAIO_MAX_SERVERS];

static flukaio_connection_t *get_connection(int cid);
static flukaio_server_t *get_server(int serverid);

/** particle buffer for ONLY one connection */
static particle_info_t  *pbuffer=NULL;
static int		 pbuffer_pos=0;
static int		 pbuffer_max=0;

int store_connection(flukaio_connection_t *conn);
int store_server(flukaio_server_t *server);

int ntrecv_(
		const int *cid,
		uint8_t *type,
		uint32_t *id, uint32_t *gen,
		double *wgt,
		double  *x, double  *y, double *z,
		double *tx, double *ty, double *tz,
		uint16_t *aa, uint16_t *zz,
		double *m, double *pc,
		double *t)
{
	ssize_t n = 0;
	particle_info_t *particle = NULL;
	flukaio_message_t msg;
	flukaio_connection_t * conn = get_connection(*cid);

	n = flukaio_receive_message(conn, &msg);

	if (n <= 0)
		return -1;

	if(type) *type = msg.type;

	switch(msg.type) {
		case N_PART:
			particle = &msg.data.particle;
			if(id)  *id  = particle->id;
			if(gen) *gen = particle->gen;
			if(wgt) *wgt = particle->weight;

			if(x) *x = particle->x;
			if(y) *y = particle->y;
			if(z) *z = particle->z;

			if(tx) *tx = particle->tx;
			if(ty) *ty = particle->ty;
			if(tz) *tz = particle->tz;

			if(aa) *aa = particle->aa;
			if(zz) *zz = particle->zz;

			if(m)  *m  = particle->m;
			if(pc) *pc = particle->pc;
			if(t)  *t  = particle->t;
			break;
		case N_IPT:
			if(gen) *gen = msg.data.ipt.turn;
			if(id)  *id  = msg.data.ipt.ipt;
			break;
		case N_EOB:
		case N_END:
		case N_ERR:
			break;
		default:
			return -1;
	}

	return n;
} // ntrecv_

int ntsendp_(	const int     *cid,
		const uint32_t *id, const uint32_t *gen,
		const double  *wgt,
		const double    *x, const double    *y, const double  *z,
		const double   *tx, const double   *ty, const double *tz,
		const uint16_t *aa, const uint16_t *zz,
		const double    *m, const double   *pc,
		const double    *t)
{
	flukaio_connection_t * conn = get_connection(*cid);

	particle_info_t part;
	part.id     = *id;
	part.gen    = *gen;
	part.weight = *wgt;
	part.x      = *x;
	part.y      = *y;
	part.z      = *z;
	part.tx     = *tx;
	part.ty     = *ty;
	part.tz     = *tz;
	part.aa     = *aa;
	part.zz     = *zz;
	part.m      = *m;
	part.pc     = *pc;
	part.t      = *t;

	return flukaio_send_particle(conn, &part);
} // ntsendp_

int ntsendpbuf_(const int     *cid,
		const uint32_t *id, const uint32_t *gen,
		const double  *wgt,
		const double    *x, const double    *y, const double  *z,
		const double   *tx, const double   *ty, const double *tz,
		const uint16_t *aa, const uint16_t *zz,
		const double    *m, const double   *pc,
		const double    *t)
{
	if (pbuffer_pos >= pbuffer_max) {
		if (pbuffer_max==0) {
			pbuffer_max = 64;	// initial size
			pbuffer = (particle_info_t*)malloc(sizeof(particle_info_t)*pbuffer_max);
		} else {
			pbuffer_max <<= 1;	// double the buffer
			pbuffer = (particle_info_t*)realloc(pbuffer, sizeof(particle_info_t)*pbuffer_max);
		}
	}
	//flukaio_connection_t * conn = get_connection(*cid);

	pbuffer[pbuffer_pos].id     = *id;
	pbuffer[pbuffer_pos].gen    = *gen;
	pbuffer[pbuffer_pos].weight = *wgt;
	pbuffer[pbuffer_pos].x      = *x;
	pbuffer[pbuffer_pos].y      = *y;
	pbuffer[pbuffer_pos].z      = *z;
	pbuffer[pbuffer_pos].tx     = *tx;
	pbuffer[pbuffer_pos].ty     = *ty;
	pbuffer[pbuffer_pos].tz     = *tz;
	pbuffer[pbuffer_pos].aa     = *aa;
	pbuffer[pbuffer_pos].zz     = *zz;
	pbuffer[pbuffer_pos].m      = *m;
	pbuffer[pbuffer_pos].pc     = *pc;
	pbuffer[pbuffer_pos].t      = *t;

	return ++pbuffer_pos;
//	return flukaio_send_particle(conn, &part);
} // ntsendpbuf_

int ntsendeob_(const int *cid)
{
	flukaio_connection_t* conn = get_connection(*cid);

	if (pbuffer_pos>0) {
		int i=0;
		for (i=0; i<pbuffer_pos; i++) {
			int err = flukaio_send_particle(conn, &pbuffer[i]);
			if (err<0) return err;
		}
		pbuffer_pos = 0;
	}

	return flukaio_send_eob(conn);
} // ntsendeob_

int ntsendipt_(const int *cid, const uint32_t *turn, const uint16_t *ipt)
{
	flukaio_connection_t* conn = get_connection(*cid);
	return flukaio_send_ipt(conn, *turn, *ipt);
} // ntsendipt_

int ntsendeoc_(const int *cid)
{
	flukaio_connection_t* conn = get_connection(*cid);
	return flukaio_send_eoc(conn);
} // ntsendeoc_

int ntsendnpart_(const int *cid, const uint32_t *npart)
{
	flukaio_connection_t* conn = get_connection(*cid);
	return flukaio_send_npart(conn, *npart);
} // ntsendnpart_

int ntnpart_(const int *cid)
{
	flukaio_message_t msg;
	flukaio_connection_t* conn = get_connection(*cid);

	if (flukaio_wait_message(conn, &msg)<0)
		return -2;

	if (msg.type == N_INT && msg.data.varint.id==N_NPART) {
		return (int)msg.data.varint.value;
	} else
		return -1;
} // ntnpart_

int ntsendbrhono_(const int *cid, const double *brhono)
{
	flukaio_connection_t* conn = get_connection(*cid);
	return flukaio_send_double(conn, N_BRHONO, *brhono);
} // ntsendbrhono_

int ntbrho_(const int *cid, double *brhono)
{
	flukaio_message_t msg;
	flukaio_connection_t* conn = get_connection(*cid);

	if (flukaio_wait_message(conn, &msg)<0)
		return -2;

	if (msg.type == N_DBLE && msg.data.vardouble.id==N_BRHONO) {
		*brhono = msg.data.vardouble.value;
		return 0;
	} else {
		*brhono = 0.0;
		return -1;
	}
} // ntbrho_

int ntwait_(
		const int *cid,
		uint8_t *type,
		uint32_t *id, uint32_t *gen,
		double *wgt,
		double  *x, double  *y, double *z,
		double *tx, double *ty, double *tz,
		uint16_t *aa, uint16_t *zz,
		double *m, double *pc,
		double *t)
{
	ssize_t n = 0;
	particle_info_t *particle = NULL;
	flukaio_message_t msg;
	flukaio_connection_t * conn = get_connection(*cid);

	n = flukaio_wait_message(conn, &msg);

	if (n <= 0)
		return -1;

	if (type) *type = msg.type;

	switch (msg.type) {
		case N_PART:
			particle = &msg.data.particle;
			if(id)  *id  = particle->id;
			if(gen) *gen = particle->gen;
			if(wgt) *wgt = particle->weight;

			if(x) *x = particle->x;
			if(y) *y = particle->y;
			if(z) *z = particle->z;

			if(tx) *tx = particle->tx;
			if(ty) *ty = particle->ty;
			if(tz) *tz = particle->tz;

			if(aa) *aa = particle->aa;
			if(zz) *zz = particle->zz;

			if(m)  *m  = particle->m;
			if(pc) *pc = particle->pc;
			if(t)  *t  = particle->t;
			break;
		case N_IPT:
			if(gen) *gen = msg.data.ipt.turn;
			if(id)  *id  = msg.data.ipt.ipt;
			break;

		case N_EOB:
		case N_END:
		case N_ERR:
		case N_INT:
	        case N_DBLE:
			break;
		default:
			return -1;
	}

	return n;
} // ntwait_

int ntserver_() {
	flukaio_server_t *server = flukaio_server_create();
	if (server)
		return store_server(server);
	return -1;
} // ntserver_

int ntstart_(const int *serverid, int *port) {
	flukaio_server_t *server = get_server(*serverid);
	if(server)
		return flukaio_server_start(server, *port);
	return -1;
} // ntstart_

int ntshdwn_(const int *serverid) {
	flukaio_server_t *server = get_server(*serverid);
	if(server) {
		int ret = flukaio_server_shutdown(server);
		// Update servers array
		servers[*serverid] = NULL;
		if (pbuffer) {
			free(pbuffer);
			pbuffer_max = pbuffer_pos = 0;
		}
		return ret;
	}
	return -1;
} // ntshdwn_

int ntgetport_(const int *serverid) {
	flukaio_server_t *server = get_server(*serverid);
	if(server)
		return flukaio_server_get_port(server);
	return -1;
} // ntgetport_

int ntaccept_(const int *serverid)
{
	flukaio_server_t *server = get_server(*serverid);
	if(server) {
		int cid;
		flukaio_connection_t *conn;
		conn = flukaio_server_accept(server);

		cid = store_connection(conn);

		// No space left
		if(cid < 0)
			connection_destroy(conn);

		return cid;
	}
	return -1;
} // ntaccept_

int ntend_(const int *cid) {
	flukaio_connection_t * conn;

	conn = get_connection(*cid);

	if (conn == NULL)
		return -1;

	flukaio_disconnect(conn);

	// Update connections array
	connections[*cid] = NULL;

	return 0;
} // ntend_

int ntrtimeout_(const int *cid, const int *seconds) {
	flukaio_connection_t * conn;

	conn = get_connection(*cid);

	if (conn == NULL)
		return -1;

	return connection_set_read_timeout(conn, *seconds);
} // ntrtimeout_

int ntwtimeout_(const int *cid, const int *seconds) {
	flukaio_connection_t * conn;

	conn = get_connection(*cid);

	if (conn == NULL)
		return -1;

	return connection_set_write_timeout(conn, *seconds);
} // ntwtimeout_

int ntconnect_(char *host, int *port, const long hostlen)
{
	int cid;
	char *chost;
	chost = create_fortran_string(host, hostlen);

	if (chost == NULL)
		return -1;

	flukaio_connection_t *conn = flukaio_connect(connection_create(-1), chost, *port);
	if (conn == NULL)
		return -1;

	free(chost);

	if (conn == NULL)
		return -1;

	cid = store_connection(conn);
	if (cid < 0) // No space left, destroy connection
		connection_destroy(conn);

	return cid;
} // ntconnect_

/**
 * Initializes the list of connections and servers
 */
int ntinit_()
{
	int i;
	for (i = 0; i < FLUKAIO_MAX_CONNECTIONS; i++) {
		connections[i] = NULL;
	}

	for (i = 0; i < FLUKAIO_MAX_SERVERS; i++) {
		servers[i] = NULL;
	}

	return 0;
} // ntinit_

/**
 * Convert a fortran string to c string
 */
char* create_fortran_string(const char *str, const long str_len) {
	long len;
	char *rslt = NULL;

	for(len = str_len-1; len >= 0; len--) {
		if (str[len] != ' ')
			break;
	}
	len++;

	if (len < 0)
		return NULL;

	rslt = malloc(len+1);
	strncpy(rslt, str, len);
	rslt[len] = 0;

	return rslt;
} // create_fortran_string

/**
 * Find connection data associated with cid
 * @return connection object
 */
static flukaio_connection_t *get_connection(int cid)
{
	return connections[cid];
} // get_connection

/**
 * Stores connection information an assigns it a connection id
 * @return connection id, -1 if error
 */
int store_connection(flukaio_connection_t *conn)
{
	int i;

	if (conn == NULL)
		return -1;

	// Find an empty slot in connections to store connection data
	for (i = 0; i < FLUKAIO_MAX_CONNECTIONS; i++) {
		if (connections[i] == NULL)
		{
			connections[i] = conn;
			return i;
		}
	}

	return -1;
} // store_connection

/**
 * Find server data associated with serverid
 * @return connection object
 */
static flukaio_server_t *get_server(int serverid)
{
	return servers[serverid];
} // get_server

/**
 * Stores server information an assigns it a server id
 * @return server id, -1 if error
 */
int store_server(flukaio_server_t *server)
{
	int i;

	if (server == NULL)
		return -1;

	// Find an empty slot in servers to store server data
	for (i = 0; i < FLUKAIO_MAX_SERVERS; i++) {
		if (servers[i] == NULL)
		{
			servers[i] = server;
			return i;
		}

	}

	return -1;
} // store_server
