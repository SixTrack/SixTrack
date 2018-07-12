#ifndef FLUKAIO_NETIO_PRIV_H__
#define FLUKAIO_NETIO_PRIV_H__

#include <stdlib.h>
#include <stdio.h>
#include <unistd.h>
#include <errno.h>
#include <string.h>
#include <sys/types.h>
#include <sys/socket.h>
#include <arpa/inet.h>
#include <fcntl.h>
#include <netinet/in.h>
#include <netinet/tcp.h>
#include <netdb.h>

ssize_t netio_read(int fd, void *buf, size_t n);
ssize_t netio_write(int fd, const void *buf, size_t len);
ssize_t netio_write_str(int fd, const char *msg);
int     netio_can_write(int fd, long timeout_sec, long timeout_usec);
int     netio_can_read(int fd, long timeout_sec, long timeout_usec);
int     netio_connect(const char *host, int portnum);
int     netio_set_nonblocking(int fd);
int     netio_set_nodelay(int fd);

int     netio_server_start(int *port, int *sockfd);
int     netio_server_accept(int server_sockfd);
int     netio_server_shutdown(int server_sockfd);

#endif
