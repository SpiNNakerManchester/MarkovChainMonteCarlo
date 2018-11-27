#include <spin1_api.h>
#include <stdint.h>
#include <stdbool.h>
#include <debug.h>
#include <data_specification.h>

#ifndef use
#define use(x) do {} while ((x)!=(x))
#endif

// Timeout between sending and receiving results in number of timer ticks
#define TIMEOUT 3

// The parameters to be read from memory
enum params {
    DATA_SIZE = 0,
    N_CHIPS,
    KEY,
    WINDOW_SIZE,
    SEQUENCE_MASK,
    TIMER,
    DATA
};

// An array of sequence numbers received on each core
// Note that this potentially will be in SDRAM with enough cores
uint *sequence_received;

// The list of keys for each of the cores - ordered for quick searching
uint *chip_keys;

// The size of the remaining data to be sent
uint data_size;

// The number of chips the data is to be sent to
uint n_chips;

// The base key to send the data with
uint key;

// The window size for sending the data
uint window_size;

// The mask which indicates the sequence number
uint sequence_mask;

// Pointer to the start of the data still to be sent and confirmed
uint *data;

// The timer tick at which the send will have timed out
uint send_timeout;

// The next sequence number to be sent
uint next_sequence = 0;

// The end sequence ignoring the wrap around of sequences
uint next_end_sequence_unwrapped = 0xFFFFFFFF;

// The end sequence with wrap around of sequences
uint next_end_sequence = 0;


void send_callback(uint send_time, uint unused) {
    use(unused);

    // If the data has all been sent, send a start message and quit
    if (data_size == 0) {
        // log_info("All data has been sent and confirmed");
        io_printf(IO_BUF, "All data has been sent and confirmed");
        while (!spin1_send_mc_packet(key, 0, NO_PAYLOAD)) {
            spin1_delay_us(1);
        }
        spin1_exit(0);
        return;
    }

    // Send the next packets
    uint packets_to_send = window_size;
    if (window_size > data_size) {
        packets_to_send = data_size;
    }
    for (uint i = 0; i < packets_to_send; i++) {
        uint sequence = (next_sequence + i) & sequence_mask;
        while (!spin1_send_mc_packet(
                key + sequence, data[i], WITH_PAYLOAD)) {
            spin1_delay_us(1);
        }
    }
    next_end_sequence_unwrapped = next_sequence + packets_to_send - 1;
    next_end_sequence = next_end_sequence_unwrapped & sequence_mask;
    send_timeout = send_time + TIMEOUT;
}


void timer_callback(uint time, uint unused) {
    use(unused);

    // Only do anything if there is sequence to check
    if (next_end_sequence_unwrapped == 0xFFFFFFFF) {
        return;
    }

    // Determine if the timeout would expire
    uint timed_out = 0;
    if (time >= send_timeout) {
        timed_out = 1;
    }

    // Find the smallest sequence number received
    uint smallest_sequence = 0xFFFFFFFF;
    for (uint i = 0; i < n_chips; i++) {
        if (sequence_received[i] < smallest_sequence) {
            smallest_sequence = sequence_received[i];

            // If one of the cores hasn't yet got the latest sequence
            // and we haven't timed out, we may as well give up
            if ((smallest_sequence < next_end_sequence_unwrapped) &&
                    !timed_out) {
                break;
            }
        }
    }

    // If all chips have the data, or we have timed out, start sending again
    if ((smallest_sequence == next_end_sequence_unwrapped) || timed_out) {

        // Adjust the sequence to the next sequence to send
        if (smallest_sequence >= next_sequence &&
                smallest_sequence <= next_end_sequence_unwrapped) {

            // Work out how many words have been sent
            uint words_sent = (smallest_sequence - next_sequence) + 1;

            // Adjust the pointers to the next bit of data to send
            data = &(data[words_sent]);
            data_size -= words_sent;

            next_sequence = (smallest_sequence + 1) & sequence_mask;
        }

        // Avoid sending more data in case the send takes too long
        next_end_sequence_unwrapped = 0xFFFFFFFF;

        // Send the data
        spin1_schedule_callback(send_callback, time, 0, 1);
    }
}

void multicast_callback(uint key, uint payload) {

    // Find the core that this packet is from using binary search
    uint imin = 0;
    uint imax = n_chips;
    while (imin < imax) {

        uint imid = (imax + imin) >> 1;

        // If the key is found, update the sequence with the payload
        if (chip_keys[imid] == key) {

            uint sequence = payload;
            uint current_sequence = sequence_received[imid];

            // If the range is wrapped, adjust the sequence to ignore the wrap
            if (next_end_sequence < next_sequence) {
                if (sequence < next_sequence) {
                    sequence += sequence_mask + 1;
                }
            }


            // Only update if the sequence is in the expected range
            if (sequence >= next_sequence &&
                    sequence <= next_end_sequence_unwrapped) {

                // If the current sequence is out of range, update
                if (current_sequence < next_sequence ||
                        current_sequence > next_end_sequence_unwrapped ||
                        sequence > current_sequence) {
                    sequence_received[imid] = sequence;
                }
            }
            break;
        }

        if (chip_keys[imid] < key) {
            imin = imid + 1;
        } else {
            imax = imid;
        }
    }
}


void empty_multicast_callback(uint key, uint payload) {
    use(key);
    use(payload);
}

void c_main() {

    address_t data_address = data_specification_get_data_address();
    address_t params = data_specification_get_region(0, data_address);

    // Get the size of the data in words
    data_size = params[DATA_SIZE];
    // log_info("Data size = %d", data_size);
    io_printf(IO_BUF, "Data size = %d", data_size);

    // Get a count of the chips to load on to
    n_chips = params[N_CHIPS];
    // log_info("N chips = %d", n_chips);
    io_printf(IO_BUF, "N chips = %d", n_chips);

    // Get the key to send the data with
    key = params[KEY];
    // log_info("Key = 0x%08x", key);
    io_printf(IO_BUF, "Key = 0x%08x", key);

    // Get the number of packets to be sent at the same time
    window_size = params[WINDOW_SIZE];
    // log_info("Window size = %d", window_size);
    io_printf(IO_BUF, "Window size = %d", window_size);

    // Get the total number of window spaces available
    sequence_mask = params[SEQUENCE_MASK];
    // log_info("Sequence mask = 0x%08x", sequence_mask);
    io_printf(IO_BUF, "Sequence mask = 0x%08x", sequence_mask);

    // Get the timer tick
    uint timer = params[TIMER];
    // log_info("Timer = %d", timer);
    io_printf(IO_BUF, "Timer = %d", timer);

    // Get a pointer to the data - not worth copying at present
    data = (uint *) &(params[DATA]);

    // Get a pointer to the keys for each chip which will verify
    // the reception of data
    uint *chip_data = (uint *) &(data[data_size]);

    // Try to allocate an array of keys in DTCM
    chip_keys = (uint *) spin1_malloc(n_chips * sizeof(uint));
    if (chip_keys == NULL) {
        log_warning("Could not allocate chip keys in DTCM - using SDRAM");
        chip_keys = chip_data;
    } else {
        spin1_memcpy(chip_keys, chip_data, n_chips * sizeof(uint));
    }

    // Try to allocate an array of last sequences in DTCM
    sequence_received = (uint *) spin1_malloc(n_chips * sizeof(uint));
    if (sequence_received == NULL) {
        log_warning("Could not allocate sequences in DTCM - using SDRAM");
        sequence_received = (uint *) sark_xalloc(
                sv->sdram_heap, n_chips * sizeof(uint), 0, ALLOC_LOCK);
        if (sequence_received == NULL) {
            log_error("Could not allocate sequences in SDRAM");
            rt_error(RTE_SWERR);
        }
    }
    for (uint i = 0; i < n_chips; i++) {
        sequence_received[i] = sequence_mask;
    }
    send_timeout = TIMEOUT;

    // Set up the timer
    spin1_set_timer_tick(timer);

    // Setup callback on multicast packet with payload, for acknowledge packets
    spin1_callback_on(MCPL_PACKET_RECEIVED, multicast_callback, -1);
    spin1_callback_on(MC_PACKET_RECEIVED, empty_multicast_callback, -1);

    // Setup callback on timer to timeout sent packets and send next packets
    spin1_callback_on(TIMER_TICK, timer_callback, 2);

    // Schedule the first call of the send callback
    spin1_schedule_callback(send_callback, 0, 0, 1);

    // Start in sync with all the cores, to ensure they are ready
    // to receive packets
    spin1_start(SYNC_WAIT);

}
