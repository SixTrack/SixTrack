#include "NetIO.h"

static ssize_t netio_xwrite(int fd, const void *buf, size_t len);

/**
 * Read if possible
 * @param fd file/socket descriptor
 * @param buf buffer
 * @param len buffer length
 * @return bytes read, -1 on failure
 */
ssize_t netio_read(const int fd, void *buf, size_t len)
{
	ssize_t nr;
	while (1) {
		nr = read(fd, buf, len);
		//if ((nr < 0) && (errno == EAGAIN || errno == EINTR))
		if ((nr < 0) && (errno == EINTR))
			continue;
		return nr;
	}
}

/**
 * Writes the complete buffer safely
 * Blocks until done
 *
 * @param fd file/socket descriptor
 * @param buf buffer
 * @param count number of bytes to write
 * @return bytes written, -1 on failure
 */
ssize_t netio_write(int fd, const void *buf, size_t count)
{
	const char *p = buf;
	ssize_t total = 0;

	while (count > 0) {
		ssize_t written = netio_xwrite(fd, p, count);
		if (written < 0)
			return -1;
		if (!written) {
			errno = ENOSPC;
			return -1;
		}
		count -= written;
		p += written;
		total += written;
	}

	return total;
}

/**
 * Writes to fd and restarts operations if recoverable error.
 *
 * Either sends something or the connection is broken
 * Doesn't guarantee len bytes are written
 *
 * Uses netio_can_write after first error happens
 * @return bytes written, -1 on failure
 */
static ssize_t netio_xwrite(int fd, const void *buf, size_t len)
{
	ssize_t nr;
	while (1) {
		nr = write(fd, buf, len);
		if ((nr < 0) && (errno == EAGAIN || errno == EINTR))
		{
			if (netio_can_write(fd, 5, 0) < 0)
				return -1;
			continue;
		}
		return nr;
	}
}


/**
 * Check if the socket can write data
 * @param fd socket descriptor to check
 * @param timeout_sec timeout in seconds
 * @param timeout_usec timeout in micro-seconds
 * @return 1 if socket available to send data
 */
int netio_can_write(int fd, long timeout_sec, long timeout_usec)
{
	int ret;
	fd_set writefds;
	struct timeval tv;

	FD_ZERO(&writefds);
	FD_SET(fd, &writefds);

	while(1) {
		tv.tv_sec  = timeout_sec;
		tv.tv_usec = timeout_usec;

		ret = select(fd+1, NULL, &writefds, NULL, &tv);
		if (ret < 0) {
			if (errno == EINTR)
				continue;
			return -1;
		}
		return FD_ISSET(fd, &writefds) != 0;
	}
}

/**
 * Check if the socket can read data
 * @param fd socket descriptor to check
 * @param timeout_sec timeout in seconds
 * @param timeout_usec timeout in micro-seconds
 * @return 1 if socket available to send data
 */
int netio_can_read(int fd, long timeout_sec, long timeout_usec)
{
	int ret;
	fd_set readfds;
	struct timeval tv;

	FD_ZERO(&readfds);
	FD_SET(fd, &readfds);

	while(1) {
		tv.tv_sec  = timeout_sec;
		tv.tv_usec = timeout_usec;

		ret = select(fd+1, &readfds, NULL, NULL, &tv);
		if (ret < 0) {
			if (errno == EINTR)
				continue;
			return -1;
		}
		return FD_ISSET(fd, &readfds) != 0;
	}
}

/**
 * Enable no delay flag in a socket
 * @return setsockopt value
 */
int netio_set_nodelay(int fd) {
	int flag = 1;
	return setsockopt(fd, IPPROTO_TCP, TCP_NODELAY, (char *)&flag, sizeof(int));
}

/**
 * Connect to host and port
 * @return 1 socket fd if ok, -1 if failure
 */
int netio_connect(const char *host, int portnum)
{
	int sockfd = -1;
	struct sockaddr_in sin;
	struct hostent *server;

	sockfd = socket(PF_INET, SOCK_STREAM, IPPROTO_TCP);
	if (sockfd < 0)
	{
		fprintf(stderr, "Unable to create socket (%s)\n", strerror(errno));
		return -1;
	}

	// Set connection parameters
	server = gethostbyname(host);
	if (server == NULL)
	{
		fprintf(stderr, "Unable to resolve host %s\n", host);
		close(sockfd);
		sockfd = -1;
		return -1;
	}

	memset(&sin, 0, sizeof(sin));
	sin.sin_family = AF_INET;
	sin.sin_port = htons(portnum);
	memcpy(&sin.sin_addr, server->h_addr, server->h_length);

	if (connect(sockfd, (struct sockaddr *)&sin, sizeof(sin)) < 0)
	{
		fprintf(stderr, "Unable to connect to %s %d error=(%s)\n", host, portnum, strerror(errno));
		close(sockfd);
		sockfd = -1;
		return -1;
	}
	netio_set_nonblocking(sockfd);
#ifdef FLUKAIO_NODELAY
	netio_set_nodelay(sockfd);
#endif

	return sockfd;
}

/**
 * Set non blocking flag in a socket
 * @return fnctl return value
 */
int netio_set_nonblocking(int fd) {
	int flags;
	if (-1 == (flags = fcntl(fd, F_GETFL, 0)))
		flags = 0;
	return fcntl(fd, F_SETFL, flags | O_NONBLOCK);
}

/**
 * Starts server and returns socket
 * @param port number where to listen (random if 0)
 * @return socket or -1 if failure
 */
int netio_server_start(int *port, int *sockfd)
{
	struct sockaddr_in myaddr_in;
	socklen_t myaddr_in_len = sizeof(myaddr_in);

	/* Create the listen socket. */
	*sockfd = socket(AF_INET, SOCK_STREAM, IPPROTO_TCP);
	if (*sockfd < 0) {
		perror("netio");
		return -1;
	}

	/* clear out address structures */
	memset(&myaddr_in, 0, sizeof(myaddr_in));

	myaddr_in.sin_family = AF_INET;
	myaddr_in.sin_port = htons(*port);	// If the user supplied 0 then it will be assigned by the kernel
	myaddr_in.sin_addr.s_addr = htonl(INADDR_ANY);

	/* Bind the listen address to the socket. */
	if (bind(*sockfd, (const struct sockaddr*)&myaddr_in, myaddr_in_len) < 0) {
		close(*sockfd);
		perror("netio");
		return -1;
	}

	/* Find the assigned port and return it back */
	if (getsockname(*sockfd, (struct sockaddr *)&myaddr_in, &myaddr_in_len) < 0) {
		close(*sockfd);
		perror("netio");
		return -1;
	}

	*port = ntohs(myaddr_in.sin_port);

	if (listen(*sockfd, 5) == -1) {
		close(*sockfd);
		perror("netio");
		return -1;
	}

	return *port;
}

/**
 * Accepts a connection in a server socket
 * @param server_sockfd an already initialized server socket
 * @return new connection socket, -1 if failure
 */
int netio_server_accept(int server_sockfd)
{
	int sockfd;
	struct sockaddr_storage peeraddr_in;
	socklen_t addrlen = sizeof(peeraddr_in);

	if (server_sockfd < 0)
		return -1;

	memset((char *)&peeraddr_in, 0, sizeof(peeraddr_in));

	sockfd = accept(server_sockfd,
			(struct sockaddr*)&peeraddr_in,
			&addrlen);

	if (sockfd < 0)
		return -1;

#ifdef FLUKAIO_NODELAY
	netio_set_nodelay(sockfd);
#endif
	netio_set_nonblocking(sockfd);

	return sockfd;
}

/**
 * Shutdownd and closes server
 * @param server_sockfd an already initialized server socket
 * @return 0 if success, -1 if failure
 */
int netio_server_shutdown(int server_sockfd)
{
	if (server_sockfd < 0)
		return 0;

	if (shutdown(server_sockfd, SHUT_RDWR) == -1) {
		perror("netio");
		return -1;
	}

	close(server_sockfd);

	return 0;
}
