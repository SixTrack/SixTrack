#include "FakeFlukaIO.h"

ssize_t fakeflukaio_insert_message(
        flukaio_connection_t *conn,
        const unsigned char type,
        const void *data, const size_t datalen)
{
    ssize_t size;
    size = connection_write_message_bin(type, data, datalen,
            conn->in_buffer + conn->in_buffer_end,
            conn->in_buffer_size - conn->in_buffer_end);

    if (size > 0) {
        conn->in_buffer_end += size;
        conn->in_buffer_len += size;
	}

    return size;
}

ssize_t fakeflukaio_insert_eob_message(flukaio_connection_t *conn) {
    return fakeflukaio_insert_message(conn, N_EOB, NULL, 0);
}

ssize_t fakeflukaio_insert_ipt_message(flukaio_connection_t *conn, const uint32_t turn, const uint16_t ipt) {
	flukaio_ipt_data_t data;
	data.turn = turn;
	data.ipt = ipt;
    return fakeflukaio_insert_message(conn, N_IPT, &data, sizeof(data));
}

ssize_t fakeflukaio_insert_particle_message(flukaio_connection_t *conn, const particle_info_t *part) {
    return fakeflukaio_insert_message(conn, N_PART, part, sizeof(particle_info_t));
}

ssize_t fakeflukaio_insert_eoc_message(flukaio_connection_t *conn) {
    return fakeflukaio_insert_message(conn, N_END, NULL, 0);
}

ssize_t fakeflukaio_insert_config_message(flukaio_connection_t *conn) {
    return fakeflukaio_insert_message(conn, N_CONF, NULL, 0);
}

ssize_t fakeflukaio_insert_unknown_message(flukaio_connection_t *conn) {
    return fakeflukaio_insert_message(conn, N_OTHER, NULL, 0);
}

ssize_t fakeflukaio_insert_hsk_message(flukaio_connection_t *conn,
		const uint16_t major, const uint16_t minor, const uint32_t key)
{
	flukaio_hsk_data_t data;
	data.minor = minor;
	data.major = major;
	data.key = key;
    return fakeflukaio_insert_message(conn, N_HSK, &data, sizeof(data));
}

ssize_t fakeflukaio_get_last_sent_message(flukaio_connection_t *conn, flukaio_message_t *msg)
{
    ssize_t len  = connection_read_message_from(msg, conn->out_buffer, conn->out_buffer_len);

    if (len > 0)
    {
        conn->out_buffer_len -= len;
        memmove(conn->out_buffer, conn->out_buffer + len, conn->out_buffer_len);
    }
    return len;
}

void fakeflukaio_clean_buffers(flukaio_connection_t *conn)
{
    conn->in_buffer_len  = 0;
    conn->out_buffer_len = 0;
}
