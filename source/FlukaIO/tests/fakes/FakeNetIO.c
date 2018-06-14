#include "FakeNetIO.h"

/* FakeNetIO */
#define DEFAULT_RAND_PORT 44
static ssize_t fakenetio_read_retval = 0;
static ssize_t fakenetio_last_sent_size = 0;
static ssize_t fakenetio_next_port = DEFAULT_RAND_PORT;


/* Function substitutes */
int fakenetio_server_start(int *port, int *sockfd) {
    if (*port == 0)
        *port = fakenetio_next_port;

    *sockfd = 101;

    fakenetio_next_port = DEFAULT_RAND_PORT;

    return *port;
}

int fakenetio_server_accept(int server_sockfd)
{
    return 102;
}

int fakenetio_server_shutdown(int server_sockfd)
{
    return 0;
}

int fakenetio_connect(const char *host, int portnum)
{
    return 103;
}

ssize_t fakenetio_write(int fd, const void *buf, size_t n)
{
    fakenetio_last_sent_size = n;
    return 0;
}

ssize_t fakenetio_read(int fd, void *buf, size_t n) {
    return fakenetio_read_retval;
}

int fakenetio_can_write(int fd, long timeout_sec, long timeout_usec)
{
    return 1;
}

int fakenetio_can_read(int fd, long timeout_sec, long timeout_usec)
{
    return 1;
}

int fakenetio_set_nonblocking(int fd) {
    return 0;
}
int fakenetio_set_nodelay(int fd) {
    return 0;
}


/* Information retreival */
ssize_t fakenetio_get_last_sent_size() {
    ssize_t n = fakenetio_last_sent_size;
    fakenetio_last_sent_size = 0;
    return n;
}

void fakenetio_set_read_retval(ssize_t n)
{
    fakenetio_read_retval = n;
}

void fakenetio_set_next_port(int port)
{
    fakenetio_next_port = port;
}
