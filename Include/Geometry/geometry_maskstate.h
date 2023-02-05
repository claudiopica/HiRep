/**
 * @file geometry_maskstate.h
 * @brief Mask structs that can give information on whether neighbors of a given site are located
 *        in the buffer and require communications first or not.
 */

#ifndef GEOMETRY_MASKSTATE_H
#define GEOMETRY_MASKSTATE_H
#ifdef __cplusplus
   extern "C" {
#endif

// this MUST fit in a char
enum MaskState {
    T_UP_MASK = (1u << 0),
    T_DN_MASK = (1u << 1),
    X_UP_MASK = (1u << 2),
    X_DN_MASK = (1u << 3),
    Y_UP_MASK = (1u << 4),
    Y_DN_MASK = (1u << 5),
    Z_UP_MASK = (1u << 6),
    Z_DN_MASK = (1u << 7),
    FULL_MASK = (1u << 8)-1
};

static inline char invertMask(char mask) {
    return mask ^ FULL_MASK;
}

#define _PRINT_BYTE "%c%c%c%c%c%c%c%c"
#define _BINARY(byte)  \
  (byte & 0x80 ? '1' : '0'), \
  (byte & 0x40 ? '1' : '0'), \
  (byte & 0x20 ? '1' : '0'), \
  (byte & 0x10 ? '1' : '0'), \
  (byte & 0x08 ? '1' : '0'), \
  (byte & 0x04 ? '1' : '0'), \
  (byte & 0x02 ? '1' : '0'), \
  (byte & 0x01 ? '1' : '0') 


#ifdef __cplusplus
   }
#endif
#endif

