#ifndef FLUKAIO_FLUKAIOSERVER_H__
#define FLUKAIO_FLUKAIOSERVER_H__

#include "Connection.h"

#ifdef __cplusplus
extern "C" {
#endif

typedef struct {
	int port;
	int fd;
	int (*start)(int *port, int *sockfd);
	int (*accept)(int server_sockfd);
	int (*shutdown)(int server_sockfd);
	flukaio_connection_t* (*connection_create)(int);
} flukaio_server_t;

flukaio_server_t* flukaio_server_create();
int flukaio_server_start(flukaio_server_t* server, unsigned int port);
flukaio_connection_t *flukaio_server_accept(flukaio_server_t* server);
int flukaio_server_shutdown(flukaio_server_t* server);

int flukaio_server_get_port(flukaio_server_t* server);

extern int (*flukaio_handshake_server)(flukaio_connection_t * conn);

#ifdef __cplusplus
}
#endif

#endif

