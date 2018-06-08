#ifndef FLUKAIO_FAKENETIO_H__
#define FLUKAIO_FAKENETIO_H__

#include <stdlib.h>

#ifdef __cplusplus
extern "C" {
#endif

/* Function substitutes */
int     fakenetio_server_start(int *port, int *sockfd);
int     fakenetio_server_accept(int server_sockfd);
int     fakenetio_server_shutdown(int server_sockfd);

int     fakenetio_connect(const char *host, int portnum);

ssize_t fakenetio_read(int fd, void *buf, size_t n);
ssize_t fakenetio_write(int fd, const void *buf, size_t len);

int     fakenetio_can_write(int fd, long timeout_sec, long timeout_usec);
int     fakenetio_can_read(int fd, long timeout_sec, long timeout_usec);

int     fakenetio_set_nonblocking(int fd);
int     fakenetio_set_nodelay(int fd);

/* Information retreival */
ssize_t fakenetio_get_last_sent_size();
void fakenetio_set_read_retval(ssize_t n);
void fakenetio_set_next_port(int port);

#ifdef __cplusplus
}
#endif

#endif
