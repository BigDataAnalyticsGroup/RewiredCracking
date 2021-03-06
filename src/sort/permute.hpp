#pragma once

#include <cstdint>
#include <immintrin.h>


const uint8_t permutation4_table[16][8] = {
    /* 0b0000 */ { 0, 1, 2, 3, 4, 5, 6, 7 },
    /* 0b0001 */ { 0, 1, 2, 3, 4, 5, 6, 7 },
    /* 0b0010 */ { 2, 3, 0, 1, 4, 5, 6, 7 },
    /* 0b0011 */ { 0, 1, 2, 3, 4, 5, 6, 7 },
    /* 0b0100 */ { 4, 5, 0, 1, 2, 3, 6, 7 },
    /* 0b0101 */ { 0, 1, 4, 5, 2, 3, 6, 7 },
    /* 0b0110 */ { 2, 3, 4, 5, 0, 1, 6, 7 },
    /* 0b0111 */ { 0, 1, 2, 3, 4, 5, 6, 7 },
    /* 0b1000 */ { 6, 7, 0, 1, 2, 3, 4, 5 },
    /* 0b1001 */ { 0, 1, 6, 7, 2, 3, 4, 5 },
    /* 0b1010 */ { 2, 3, 6, 7, 0, 1, 4, 5 },
    /* 0b1011 */ { 0, 1, 2, 3, 6, 7, 4, 5 },
    /* 0b1100 */ { 4, 5, 6, 7, 0, 1, 2, 3 },
    /* 0b1101 */ { 0, 1, 4, 5, 6, 7, 2, 3 },
    /* 0b1110 */ { 2, 3, 4, 5, 6, 7, 0, 1 },
    /* 0b1111 */ { 0, 1, 2, 3, 4, 5, 6, 7 },
};

const uint8_t permutation8_table[256][8] = {
    /* 0b00000000 */ { 0, 1, 2, 3, 4, 5, 6, 7 },
    /* 0b00000001 */ { 0, 1, 2, 3, 4, 5, 6, 7 },
    /* 0b00000010 */ { 1, 0, 2, 3, 4, 5, 6, 7 },
    /* 0b00000011 */ { 0, 1, 2, 3, 4, 5, 6, 7 },
    /* 0b00000100 */ { 2, 0, 1, 3, 4, 5, 6, 7 },
    /* 0b00000101 */ { 0, 2, 1, 3, 4, 5, 6, 7 },
    /* 0b00000110 */ { 1, 2, 0, 3, 4, 5, 6, 7 },
    /* 0b00000111 */ { 0, 1, 2, 3, 4, 5, 6, 7 },
    /* 0b00001000 */ { 3, 0, 1, 2, 4, 5, 6, 7 },
    /* 0b00001001 */ { 0, 3, 1, 2, 4, 5, 6, 7 },
    /* 0b00001010 */ { 1, 3, 0, 2, 4, 5, 6, 7 },
    /* 0b00001011 */ { 0, 1, 3, 2, 4, 5, 6, 7 },
    /* 0b00001100 */ { 2, 3, 0, 1, 4, 5, 6, 7 },
    /* 0b00001101 */ { 0, 2, 3, 1, 4, 5, 6, 7 },
    /* 0b00001110 */ { 1, 2, 3, 0, 4, 5, 6, 7 },
    /* 0b00001111 */ { 0, 1, 2, 3, 4, 5, 6, 7 },
    /* 0b00010000 */ { 4, 0, 1, 2, 3, 5, 6, 7 },
    /* 0b00010001 */ { 0, 4, 1, 2, 3, 5, 6, 7 },
    /* 0b00010010 */ { 1, 4, 0, 2, 3, 5, 6, 7 },
    /* 0b00010011 */ { 0, 1, 4, 2, 3, 5, 6, 7 },
    /* 0b00010100 */ { 2, 4, 0, 1, 3, 5, 6, 7 },
    /* 0b00010101 */ { 0, 2, 4, 1, 3, 5, 6, 7 },
    /* 0b00010110 */ { 1, 2, 4, 0, 3, 5, 6, 7 },
    /* 0b00010111 */ { 0, 1, 2, 4, 3, 5, 6, 7 },
    /* 0b00011000 */ { 3, 4, 0, 1, 2, 5, 6, 7 },
    /* 0b00011001 */ { 0, 3, 4, 1, 2, 5, 6, 7 },
    /* 0b00011010 */ { 1, 3, 4, 0, 2, 5, 6, 7 },
    /* 0b00011011 */ { 0, 1, 3, 4, 2, 5, 6, 7 },
    /* 0b00011100 */ { 2, 3, 4, 0, 1, 5, 6, 7 },
    /* 0b00011101 */ { 0, 2, 3, 4, 1, 5, 6, 7 },
    /* 0b00011110 */ { 1, 2, 3, 4, 0, 5, 6, 7 },
    /* 0b00011111 */ { 0, 1, 2, 3, 4, 5, 6, 7 },
    /* 0b00100000 */ { 5, 0, 1, 2, 3, 4, 6, 7 },
    /* 0b00100001 */ { 0, 5, 1, 2, 3, 4, 6, 7 },
    /* 0b00100010 */ { 1, 5, 0, 2, 3, 4, 6, 7 },
    /* 0b00100011 */ { 0, 1, 5, 2, 3, 4, 6, 7 },
    /* 0b00100100 */ { 2, 5, 0, 1, 3, 4, 6, 7 },
    /* 0b00100101 */ { 0, 2, 5, 1, 3, 4, 6, 7 },
    /* 0b00100110 */ { 1, 2, 5, 0, 3, 4, 6, 7 },
    /* 0b00100111 */ { 0, 1, 2, 5, 3, 4, 6, 7 },
    /* 0b00101000 */ { 3, 5, 0, 1, 2, 4, 6, 7 },
    /* 0b00101001 */ { 0, 3, 5, 1, 2, 4, 6, 7 },
    /* 0b00101010 */ { 1, 3, 5, 0, 2, 4, 6, 7 },
    /* 0b00101011 */ { 0, 1, 3, 5, 2, 4, 6, 7 },
    /* 0b00101100 */ { 2, 3, 5, 0, 1, 4, 6, 7 },
    /* 0b00101101 */ { 0, 2, 3, 5, 1, 4, 6, 7 },
    /* 0b00101110 */ { 1, 2, 3, 5, 0, 4, 6, 7 },
    /* 0b00101111 */ { 0, 1, 2, 3, 5, 4, 6, 7 },
    /* 0b00110000 */ { 4, 5, 0, 1, 2, 3, 6, 7 },
    /* 0b00110001 */ { 0, 4, 5, 1, 2, 3, 6, 7 },
    /* 0b00110010 */ { 1, 4, 5, 0, 2, 3, 6, 7 },
    /* 0b00110011 */ { 0, 1, 4, 5, 2, 3, 6, 7 },
    /* 0b00110100 */ { 2, 4, 5, 0, 1, 3, 6, 7 },
    /* 0b00110101 */ { 0, 2, 4, 5, 1, 3, 6, 7 },
    /* 0b00110110 */ { 1, 2, 4, 5, 0, 3, 6, 7 },
    /* 0b00110111 */ { 0, 1, 2, 4, 5, 3, 6, 7 },
    /* 0b00111000 */ { 3, 4, 5, 0, 1, 2, 6, 7 },
    /* 0b00111001 */ { 0, 3, 4, 5, 1, 2, 6, 7 },
    /* 0b00111010 */ { 1, 3, 4, 5, 0, 2, 6, 7 },
    /* 0b00111011 */ { 0, 1, 3, 4, 5, 2, 6, 7 },
    /* 0b00111100 */ { 2, 3, 4, 5, 0, 1, 6, 7 },
    /* 0b00111101 */ { 0, 2, 3, 4, 5, 1, 6, 7 },
    /* 0b00111110 */ { 1, 2, 3, 4, 5, 0, 6, 7 },
    /* 0b00111111 */ { 0, 1, 2, 3, 4, 5, 6, 7 },
    /* 0b01000000 */ { 6, 0, 1, 2, 3, 4, 5, 7 },
    /* 0b01000001 */ { 0, 6, 1, 2, 3, 4, 5, 7 },
    /* 0b01000010 */ { 1, 6, 0, 2, 3, 4, 5, 7 },
    /* 0b01000011 */ { 0, 1, 6, 2, 3, 4, 5, 7 },
    /* 0b01000100 */ { 2, 6, 0, 1, 3, 4, 5, 7 },
    /* 0b01000101 */ { 0, 2, 6, 1, 3, 4, 5, 7 },
    /* 0b01000110 */ { 1, 2, 6, 0, 3, 4, 5, 7 },
    /* 0b01000111 */ { 0, 1, 2, 6, 3, 4, 5, 7 },
    /* 0b01001000 */ { 3, 6, 0, 1, 2, 4, 5, 7 },
    /* 0b01001001 */ { 0, 3, 6, 1, 2, 4, 5, 7 },
    /* 0b01001010 */ { 1, 3, 6, 0, 2, 4, 5, 7 },
    /* 0b01001011 */ { 0, 1, 3, 6, 2, 4, 5, 7 },
    /* 0b01001100 */ { 2, 3, 6, 0, 1, 4, 5, 7 },
    /* 0b01001101 */ { 0, 2, 3, 6, 1, 4, 5, 7 },
    /* 0b01001110 */ { 1, 2, 3, 6, 0, 4, 5, 7 },
    /* 0b01001111 */ { 0, 1, 2, 3, 6, 4, 5, 7 },
    /* 0b01010000 */ { 4, 6, 0, 1, 2, 3, 5, 7 },
    /* 0b01010001 */ { 0, 4, 6, 1, 2, 3, 5, 7 },
    /* 0b01010010 */ { 1, 4, 6, 0, 2, 3, 5, 7 },
    /* 0b01010011 */ { 0, 1, 4, 6, 2, 3, 5, 7 },
    /* 0b01010100 */ { 2, 4, 6, 0, 1, 3, 5, 7 },
    /* 0b01010101 */ { 0, 2, 4, 6, 1, 3, 5, 7 },
    /* 0b01010110 */ { 1, 2, 4, 6, 0, 3, 5, 7 },
    /* 0b01010111 */ { 0, 1, 2, 4, 6, 3, 5, 7 },
    /* 0b01011000 */ { 3, 4, 6, 0, 1, 2, 5, 7 },
    /* 0b01011001 */ { 0, 3, 4, 6, 1, 2, 5, 7 },
    /* 0b01011010 */ { 1, 3, 4, 6, 0, 2, 5, 7 },
    /* 0b01011011 */ { 0, 1, 3, 4, 6, 2, 5, 7 },
    /* 0b01011100 */ { 2, 3, 4, 6, 0, 1, 5, 7 },
    /* 0b01011101 */ { 0, 2, 3, 4, 6, 1, 5, 7 },
    /* 0b01011110 */ { 1, 2, 3, 4, 6, 0, 5, 7 },
    /* 0b01011111 */ { 0, 1, 2, 3, 4, 6, 5, 7 },
    /* 0b01100000 */ { 5, 6, 0, 1, 2, 3, 4, 7 },
    /* 0b01100001 */ { 0, 5, 6, 1, 2, 3, 4, 7 },
    /* 0b01100010 */ { 1, 5, 6, 0, 2, 3, 4, 7 },
    /* 0b01100011 */ { 0, 1, 5, 6, 2, 3, 4, 7 },
    /* 0b01100100 */ { 2, 5, 6, 0, 1, 3, 4, 7 },
    /* 0b01100101 */ { 0, 2, 5, 6, 1, 3, 4, 7 },
    /* 0b01100110 */ { 1, 2, 5, 6, 0, 3, 4, 7 },
    /* 0b01100111 */ { 0, 1, 2, 5, 6, 3, 4, 7 },
    /* 0b01101000 */ { 3, 5, 6, 0, 1, 2, 4, 7 },
    /* 0b01101001 */ { 0, 3, 5, 6, 1, 2, 4, 7 },
    /* 0b01101010 */ { 1, 3, 5, 6, 0, 2, 4, 7 },
    /* 0b01101011 */ { 0, 1, 3, 5, 6, 2, 4, 7 },
    /* 0b01101100 */ { 2, 3, 5, 6, 0, 1, 4, 7 },
    /* 0b01101101 */ { 0, 2, 3, 5, 6, 1, 4, 7 },
    /* 0b01101110 */ { 1, 2, 3, 5, 6, 0, 4, 7 },
    /* 0b01101111 */ { 0, 1, 2, 3, 5, 6, 4, 7 },
    /* 0b01110000 */ { 4, 5, 6, 0, 1, 2, 3, 7 },
    /* 0b01110001 */ { 0, 4, 5, 6, 1, 2, 3, 7 },
    /* 0b01110010 */ { 1, 4, 5, 6, 0, 2, 3, 7 },
    /* 0b01110011 */ { 0, 1, 4, 5, 6, 2, 3, 7 },
    /* 0b01110100 */ { 2, 4, 5, 6, 0, 1, 3, 7 },
    /* 0b01110101 */ { 0, 2, 4, 5, 6, 1, 3, 7 },
    /* 0b01110110 */ { 1, 2, 4, 5, 6, 0, 3, 7 },
    /* 0b01110111 */ { 0, 1, 2, 4, 5, 6, 3, 7 },
    /* 0b01111000 */ { 3, 4, 5, 6, 0, 1, 2, 7 },
    /* 0b01111001 */ { 0, 3, 4, 5, 6, 1, 2, 7 },
    /* 0b01111010 */ { 1, 3, 4, 5, 6, 0, 2, 7 },
    /* 0b01111011 */ { 0, 1, 3, 4, 5, 6, 2, 7 },
    /* 0b01111100 */ { 2, 3, 4, 5, 6, 0, 1, 7 },
    /* 0b01111101 */ { 0, 2, 3, 4, 5, 6, 1, 7 },
    /* 0b01111110 */ { 1, 2, 3, 4, 5, 6, 0, 7 },
    /* 0b01111111 */ { 0, 1, 2, 3, 4, 5, 6, 7 },
    /* 0b10000000 */ { 7, 0, 1, 2, 3, 4, 5, 6 },
    /* 0b10000001 */ { 0, 7, 1, 2, 3, 4, 5, 6 },
    /* 0b10000010 */ { 1, 7, 0, 2, 3, 4, 5, 6 },
    /* 0b10000011 */ { 0, 1, 7, 2, 3, 4, 5, 6 },
    /* 0b10000100 */ { 2, 7, 0, 1, 3, 4, 5, 6 },
    /* 0b10000101 */ { 0, 2, 7, 1, 3, 4, 5, 6 },
    /* 0b10000110 */ { 1, 2, 7, 0, 3, 4, 5, 6 },
    /* 0b10000111 */ { 0, 1, 2, 7, 3, 4, 5, 6 },
    /* 0b10001000 */ { 3, 7, 0, 1, 2, 4, 5, 6 },
    /* 0b10001001 */ { 0, 3, 7, 1, 2, 4, 5, 6 },
    /* 0b10001010 */ { 1, 3, 7, 0, 2, 4, 5, 6 },
    /* 0b10001011 */ { 0, 1, 3, 7, 2, 4, 5, 6 },
    /* 0b10001100 */ { 2, 3, 7, 0, 1, 4, 5, 6 },
    /* 0b10001101 */ { 0, 2, 3, 7, 1, 4, 5, 6 },
    /* 0b10001110 */ { 1, 2, 3, 7, 0, 4, 5, 6 },
    /* 0b10001111 */ { 0, 1, 2, 3, 7, 4, 5, 6 },
    /* 0b10010000 */ { 4, 7, 0, 1, 2, 3, 5, 6 },
    /* 0b10010001 */ { 0, 4, 7, 1, 2, 3, 5, 6 },
    /* 0b10010010 */ { 1, 4, 7, 0, 2, 3, 5, 6 },
    /* 0b10010011 */ { 0, 1, 4, 7, 2, 3, 5, 6 },
    /* 0b10010100 */ { 2, 4, 7, 0, 1, 3, 5, 6 },
    /* 0b10010101 */ { 0, 2, 4, 7, 1, 3, 5, 6 },
    /* 0b10010110 */ { 1, 2, 4, 7, 0, 3, 5, 6 },
    /* 0b10010111 */ { 0, 1, 2, 4, 7, 3, 5, 6 },
    /* 0b10011000 */ { 3, 4, 7, 0, 1, 2, 5, 6 },
    /* 0b10011001 */ { 0, 3, 4, 7, 1, 2, 5, 6 },
    /* 0b10011010 */ { 1, 3, 4, 7, 0, 2, 5, 6 },
    /* 0b10011011 */ { 0, 1, 3, 4, 7, 2, 5, 6 },
    /* 0b10011100 */ { 2, 3, 4, 7, 0, 1, 5, 6 },
    /* 0b10011101 */ { 0, 2, 3, 4, 7, 1, 5, 6 },
    /* 0b10011110 */ { 1, 2, 3, 4, 7, 0, 5, 6 },
    /* 0b10011111 */ { 0, 1, 2, 3, 4, 7, 5, 6 },
    /* 0b10100000 */ { 5, 7, 0, 1, 2, 3, 4, 6 },
    /* 0b10100001 */ { 0, 5, 7, 1, 2, 3, 4, 6 },
    /* 0b10100010 */ { 1, 5, 7, 0, 2, 3, 4, 6 },
    /* 0b10100011 */ { 0, 1, 5, 7, 2, 3, 4, 6 },
    /* 0b10100100 */ { 2, 5, 7, 0, 1, 3, 4, 6 },
    /* 0b10100101 */ { 0, 2, 5, 7, 1, 3, 4, 6 },
    /* 0b10100110 */ { 1, 2, 5, 7, 0, 3, 4, 6 },
    /* 0b10100111 */ { 0, 1, 2, 5, 7, 3, 4, 6 },
    /* 0b10101000 */ { 3, 5, 7, 0, 1, 2, 4, 6 },
    /* 0b10101001 */ { 0, 3, 5, 7, 1, 2, 4, 6 },
    /* 0b10101010 */ { 1, 3, 5, 7, 0, 2, 4, 6 },
    /* 0b10101011 */ { 0, 1, 3, 5, 7, 2, 4, 6 },
    /* 0b10101100 */ { 2, 3, 5, 7, 0, 1, 4, 6 },
    /* 0b10101101 */ { 0, 2, 3, 5, 7, 1, 4, 6 },
    /* 0b10101110 */ { 1, 2, 3, 5, 7, 0, 4, 6 },
    /* 0b10101111 */ { 0, 1, 2, 3, 5, 7, 4, 6 },
    /* 0b10110000 */ { 4, 5, 7, 0, 1, 2, 3, 6 },
    /* 0b10110001 */ { 0, 4, 5, 7, 1, 2, 3, 6 },
    /* 0b10110010 */ { 1, 4, 5, 7, 0, 2, 3, 6 },
    /* 0b10110011 */ { 0, 1, 4, 5, 7, 2, 3, 6 },
    /* 0b10110100 */ { 2, 4, 5, 7, 0, 1, 3, 6 },
    /* 0b10110101 */ { 0, 2, 4, 5, 7, 1, 3, 6 },
    /* 0b10110110 */ { 1, 2, 4, 5, 7, 0, 3, 6 },
    /* 0b10110111 */ { 0, 1, 2, 4, 5, 7, 3, 6 },
    /* 0b10111000 */ { 3, 4, 5, 7, 0, 1, 2, 6 },
    /* 0b10111001 */ { 0, 3, 4, 5, 7, 1, 2, 6 },
    /* 0b10111010 */ { 1, 3, 4, 5, 7, 0, 2, 6 },
    /* 0b10111011 */ { 0, 1, 3, 4, 5, 7, 2, 6 },
    /* 0b10111100 */ { 2, 3, 4, 5, 7, 0, 1, 6 },
    /* 0b10111101 */ { 0, 2, 3, 4, 5, 7, 1, 6 },
    /* 0b10111110 */ { 1, 2, 3, 4, 5, 7, 0, 6 },
    /* 0b10111111 */ { 0, 1, 2, 3, 4, 5, 7, 6 },
    /* 0b11000000 */ { 6, 7, 0, 1, 2, 3, 4, 5 },
    /* 0b11000001 */ { 0, 6, 7, 1, 2, 3, 4, 5 },
    /* 0b11000010 */ { 1, 6, 7, 0, 2, 3, 4, 5 },
    /* 0b11000011 */ { 0, 1, 6, 7, 2, 3, 4, 5 },
    /* 0b11000100 */ { 2, 6, 7, 0, 1, 3, 4, 5 },
    /* 0b11000101 */ { 0, 2, 6, 7, 1, 3, 4, 5 },
    /* 0b11000110 */ { 1, 2, 6, 7, 0, 3, 4, 5 },
    /* 0b11000111 */ { 0, 1, 2, 6, 7, 3, 4, 5 },
    /* 0b11001000 */ { 3, 6, 7, 0, 1, 2, 4, 5 },
    /* 0b11001001 */ { 0, 3, 6, 7, 1, 2, 4, 5 },
    /* 0b11001010 */ { 1, 3, 6, 7, 0, 2, 4, 5 },
    /* 0b11001011 */ { 0, 1, 3, 6, 7, 2, 4, 5 },
    /* 0b11001100 */ { 2, 3, 6, 7, 0, 1, 4, 5 },
    /* 0b11001101 */ { 0, 2, 3, 6, 7, 1, 4, 5 },
    /* 0b11001110 */ { 1, 2, 3, 6, 7, 0, 4, 5 },
    /* 0b11001111 */ { 0, 1, 2, 3, 6, 7, 4, 5 },
    /* 0b11010000 */ { 4, 6, 7, 0, 1, 2, 3, 5 },
    /* 0b11010001 */ { 0, 4, 6, 7, 1, 2, 3, 5 },
    /* 0b11010010 */ { 1, 4, 6, 7, 0, 2, 3, 5 },
    /* 0b11010011 */ { 0, 1, 4, 6, 7, 2, 3, 5 },
    /* 0b11010100 */ { 2, 4, 6, 7, 0, 1, 3, 5 },
    /* 0b11010101 */ { 0, 2, 4, 6, 7, 1, 3, 5 },
    /* 0b11010110 */ { 1, 2, 4, 6, 7, 0, 3, 5 },
    /* 0b11010111 */ { 0, 1, 2, 4, 6, 7, 3, 5 },
    /* 0b11011000 */ { 3, 4, 6, 7, 0, 1, 2, 5 },
    /* 0b11011001 */ { 0, 3, 4, 6, 7, 1, 2, 5 },
    /* 0b11011010 */ { 1, 3, 4, 6, 7, 0, 2, 5 },
    /* 0b11011011 */ { 0, 1, 3, 4, 6, 7, 2, 5 },
    /* 0b11011100 */ { 2, 3, 4, 6, 7, 0, 1, 5 },
    /* 0b11011101 */ { 0, 2, 3, 4, 6, 7, 1, 5 },
    /* 0b11011110 */ { 1, 2, 3, 4, 6, 7, 0, 5 },
    /* 0b11011111 */ { 0, 1, 2, 3, 4, 6, 7, 5 },
    /* 0b11100000 */ { 5, 6, 7, 0, 1, 2, 3, 4 },
    /* 0b11100001 */ { 0, 5, 6, 7, 1, 2, 3, 4 },
    /* 0b11100010 */ { 1, 5, 6, 7, 0, 2, 3, 4 },
    /* 0b11100011 */ { 0, 1, 5, 6, 7, 2, 3, 4 },
    /* 0b11100100 */ { 2, 5, 6, 7, 0, 1, 3, 4 },
    /* 0b11100101 */ { 0, 2, 5, 6, 7, 1, 3, 4 },
    /* 0b11100110 */ { 1, 2, 5, 6, 7, 0, 3, 4 },
    /* 0b11100111 */ { 0, 1, 2, 5, 6, 7, 3, 4 },
    /* 0b11101000 */ { 3, 5, 6, 7, 0, 1, 2, 4 },
    /* 0b11101001 */ { 0, 3, 5, 6, 7, 1, 2, 4 },
    /* 0b11101010 */ { 1, 3, 5, 6, 7, 0, 2, 4 },
    /* 0b11101011 */ { 0, 1, 3, 5, 6, 7, 2, 4 },
    /* 0b11101100 */ { 2, 3, 5, 6, 7, 0, 1, 4 },
    /* 0b11101101 */ { 0, 2, 3, 5, 6, 7, 1, 4 },
    /* 0b11101110 */ { 1, 2, 3, 5, 6, 7, 0, 4 },
    /* 0b11101111 */ { 0, 1, 2, 3, 5, 6, 7, 4 },
    /* 0b11110000 */ { 4, 5, 6, 7, 0, 1, 2, 3 },
    /* 0b11110001 */ { 0, 4, 5, 6, 7, 1, 2, 3 },
    /* 0b11110010 */ { 1, 4, 5, 6, 7, 0, 2, 3 },
    /* 0b11110011 */ { 0, 1, 4, 5, 6, 7, 2, 3 },
    /* 0b11110100 */ { 2, 4, 5, 6, 7, 0, 1, 3 },
    /* 0b11110101 */ { 0, 2, 4, 5, 6, 7, 1, 3 },
    /* 0b11110110 */ { 1, 2, 4, 5, 6, 7, 0, 3 },
    /* 0b11110111 */ { 0, 1, 2, 4, 5, 6, 7, 3 },
    /* 0b11111000 */ { 3, 4, 5, 6, 7, 0, 1, 2 },
    /* 0b11111001 */ { 0, 3, 4, 5, 6, 7, 1, 2 },
    /* 0b11111010 */ { 1, 3, 4, 5, 6, 7, 0, 2 },
    /* 0b11111011 */ { 0, 1, 3, 4, 5, 6, 7, 2 },
    /* 0b11111100 */ { 2, 3, 4, 5, 6, 7, 0, 1 },
    /* 0b11111101 */ { 0, 2, 3, 4, 5, 6, 7, 1 },
    /* 0b11111110 */ { 1, 2, 3, 4, 5, 6, 7, 0 },
    /* 0b11111111 */ { 0, 1, 2, 3, 4, 5, 6, 7 },
};

const uint8_t permutation4_packable_table[16][8] = {
    /* 0b0000 */ { 0, 1, 4, 5, 2, 3, 6, 7 },
    /* 0b0001 */ { 0, 1, 4, 5, 2, 3, 6, 7 },
    /* 0b0010 */ { 2, 3, 4, 5, 0, 1, 6, 7 },
    /* 0b0011 */ { 0, 1, 4, 5, 2, 3, 6, 7 },
    /* 0b0100 */ { 4, 5, 2, 3, 0, 1, 6, 7 },
    /* 0b0101 */ { 0, 1, 2, 3, 4, 5, 6, 7 },
    /* 0b0110 */ { 2, 3, 0, 1, 4, 5, 6, 7 },
    /* 0b0111 */ { 0, 1, 4, 5, 2, 3, 6, 7 },
    /* 0b1000 */ { 6, 7, 2, 3, 0, 1, 4, 5 },
    /* 0b1001 */ { 0, 1, 2, 3, 6, 7, 4, 5 },
    /* 0b1010 */ { 2, 3, 0, 1, 6, 7, 4, 5 },
    /* 0b1011 */ { 0, 1, 6, 7, 2, 3, 4, 5 },
    /* 0b1100 */ { 4, 5, 0, 1, 6, 7, 2, 3 },
    /* 0b1101 */ { 0, 1, 6, 7, 4, 5, 2, 3 },
    /* 0b1110 */ { 2, 3, 6, 7, 4, 5, 0, 1 },
    /* 0b1111 */ { 0, 1, 4, 5, 2, 3, 6, 7 },
};

const uint8_t permutation8_packable_table[256][8] = {
    /* 0b00000000 */ { 0, 1, 4, 5, 2, 3, 6, 7 },
    /* 0b00000001 */ { 0, 1, 4, 5, 2, 3, 6, 7 },
    /* 0b00000010 */ { 1, 0, 4, 5, 2, 3, 6, 7 },
    /* 0b00000011 */ { 0, 1, 4, 5, 2, 3, 6, 7 },
    /* 0b00000100 */ { 2, 0, 4, 5, 1, 3, 6, 7 },
    /* 0b00000101 */ { 0, 2, 4, 5, 1, 3, 6, 7 },
    /* 0b00000110 */ { 1, 2, 4, 5, 0, 3, 6, 7 },
    /* 0b00000111 */ { 0, 1, 4, 5, 2, 3, 6, 7 },
    /* 0b00001000 */ { 3, 0, 4, 5, 1, 2, 6, 7 },
    /* 0b00001001 */ { 0, 3, 4, 5, 1, 2, 6, 7 },
    /* 0b00001010 */ { 1, 3, 4, 5, 0, 2, 6, 7 },
    /* 0b00001011 */ { 0, 1, 4, 5, 3, 2, 6, 7 },
    /* 0b00001100 */ { 2, 3, 4, 5, 0, 1, 6, 7 },
    /* 0b00001101 */ { 0, 2, 4, 5, 3, 1, 6, 7 },
    /* 0b00001110 */ { 1, 2, 4, 5, 3, 0, 6, 7 },
    /* 0b00001111 */ { 0, 1, 4, 5, 2, 3, 6, 7 },
    /* 0b00010000 */ { 4, 0, 3, 5, 1, 2, 6, 7 },
    /* 0b00010001 */ { 0, 4, 3, 5, 1, 2, 6, 7 },
    /* 0b00010010 */ { 1, 4, 3, 5, 0, 2, 6, 7 },
    /* 0b00010011 */ { 0, 1, 3, 5, 4, 2, 6, 7 },
    /* 0b00010100 */ { 2, 4, 3, 5, 0, 1, 6, 7 },
    /* 0b00010101 */ { 0, 2, 3, 5, 4, 1, 6, 7 },
    /* 0b00010110 */ { 1, 2, 3, 5, 4, 0, 6, 7 },
    /* 0b00010111 */ { 0, 1, 3, 5, 2, 4, 6, 7 },
    /* 0b00011000 */ { 3, 4, 2, 5, 0, 1, 6, 7 },
    /* 0b00011001 */ { 0, 3, 2, 5, 4, 1, 6, 7 },
    /* 0b00011010 */ { 1, 3, 2, 5, 4, 0, 6, 7 },
    /* 0b00011011 */ { 0, 1, 2, 5, 3, 4, 6, 7 },
    /* 0b00011100 */ { 2, 3, 1, 5, 4, 0, 6, 7 },
    /* 0b00011101 */ { 0, 2, 1, 5, 3, 4, 6, 7 },
    /* 0b00011110 */ { 1, 2, 0, 5, 3, 4, 6, 7 },
    /* 0b00011111 */ { 0, 1, 4, 5, 2, 3, 6, 7 },
    /* 0b00100000 */ { 5, 0, 3, 4, 1, 2, 6, 7 },
    /* 0b00100001 */ { 0, 5, 3, 4, 1, 2, 6, 7 },
    /* 0b00100010 */ { 1, 5, 3, 4, 0, 2, 6, 7 },
    /* 0b00100011 */ { 0, 1, 3, 4, 5, 2, 6, 7 },
    /* 0b00100100 */ { 2, 5, 3, 4, 0, 1, 6, 7 },
    /* 0b00100101 */ { 0, 2, 3, 4, 5, 1, 6, 7 },
    /* 0b00100110 */ { 1, 2, 3, 4, 5, 0, 6, 7 },
    /* 0b00100111 */ { 0, 1, 3, 4, 2, 5, 6, 7 },
    /* 0b00101000 */ { 3, 5, 2, 4, 0, 1, 6, 7 },
    /* 0b00101001 */ { 0, 3, 2, 4, 5, 1, 6, 7 },
    /* 0b00101010 */ { 1, 3, 2, 4, 5, 0, 6, 7 },
    /* 0b00101011 */ { 0, 1, 2, 4, 3, 5, 6, 7 },
    /* 0b00101100 */ { 2, 3, 1, 4, 5, 0, 6, 7 },
    /* 0b00101101 */ { 0, 2, 1, 4, 3, 5, 6, 7 },
    /* 0b00101110 */ { 1, 2, 0, 4, 3, 5, 6, 7 },
    /* 0b00101111 */ { 0, 1, 5, 4, 2, 3, 6, 7 },
    /* 0b00110000 */ { 4, 5, 2, 3, 0, 1, 6, 7 },
    /* 0b00110001 */ { 0, 4, 2, 3, 5, 1, 6, 7 },
    /* 0b00110010 */ { 1, 4, 2, 3, 5, 0, 6, 7 },
    /* 0b00110011 */ { 0, 1, 2, 3, 4, 5, 6, 7 },
    /* 0b00110100 */ { 2, 4, 1, 3, 5, 0, 6, 7 },
    /* 0b00110101 */ { 0, 2, 1, 3, 4, 5, 6, 7 },
    /* 0b00110110 */ { 1, 2, 0, 3, 4, 5, 6, 7 },
    /* 0b00110111 */ { 0, 1, 5, 3, 2, 4, 6, 7 },
    /* 0b00111000 */ { 3, 4, 1, 2, 5, 0, 6, 7 },
    /* 0b00111001 */ { 0, 3, 1, 2, 4, 5, 6, 7 },
    /* 0b00111010 */ { 1, 3, 0, 2, 4, 5, 6, 7 },
    /* 0b00111011 */ { 0, 1, 5, 2, 3, 4, 6, 7 },
    /* 0b00111100 */ { 2, 3, 0, 1, 4, 5, 6, 7 },
    /* 0b00111101 */ { 0, 2, 5, 1, 3, 4, 6, 7 },
    /* 0b00111110 */ { 1, 2, 5, 0, 3, 4, 6, 7 },
    /* 0b00111111 */ { 0, 1, 4, 5, 2, 3, 6, 7 },
    /* 0b01000000 */ { 6, 0, 3, 4, 1, 2, 5, 7 },
    /* 0b01000001 */ { 0, 6, 3, 4, 1, 2, 5, 7 },
    /* 0b01000010 */ { 1, 6, 3, 4, 0, 2, 5, 7 },
    /* 0b01000011 */ { 0, 1, 3, 4, 6, 2, 5, 7 },
    /* 0b01000100 */ { 2, 6, 3, 4, 0, 1, 5, 7 },
    /* 0b01000101 */ { 0, 2, 3, 4, 6, 1, 5, 7 },
    /* 0b01000110 */ { 1, 2, 3, 4, 6, 0, 5, 7 },
    /* 0b01000111 */ { 0, 1, 3, 4, 2, 6, 5, 7 },
    /* 0b01001000 */ { 3, 6, 2, 4, 0, 1, 5, 7 },
    /* 0b01001001 */ { 0, 3, 2, 4, 6, 1, 5, 7 },
    /* 0b01001010 */ { 1, 3, 2, 4, 6, 0, 5, 7 },
    /* 0b01001011 */ { 0, 1, 2, 4, 3, 6, 5, 7 },
    /* 0b01001100 */ { 2, 3, 1, 4, 6, 0, 5, 7 },
    /* 0b01001101 */ { 0, 2, 1, 4, 3, 6, 5, 7 },
    /* 0b01001110 */ { 1, 2, 0, 4, 3, 6, 5, 7 },
    /* 0b01001111 */ { 0, 1, 6, 4, 2, 3, 5, 7 },
    /* 0b01010000 */ { 4, 6, 2, 3, 0, 1, 5, 7 },
    /* 0b01010001 */ { 0, 4, 2, 3, 6, 1, 5, 7 },
    /* 0b01010010 */ { 1, 4, 2, 3, 6, 0, 5, 7 },
    /* 0b01010011 */ { 0, 1, 2, 3, 4, 6, 5, 7 },
    /* 0b01010100 */ { 2, 4, 1, 3, 6, 0, 5, 7 },
    /* 0b01010101 */ { 0, 2, 1, 3, 4, 6, 5, 7 },
    /* 0b01010110 */ { 1, 2, 0, 3, 4, 6, 5, 7 },
    /* 0b01010111 */ { 0, 1, 6, 3, 2, 4, 5, 7 },
    /* 0b01011000 */ { 3, 4, 1, 2, 6, 0, 5, 7 },
    /* 0b01011001 */ { 0, 3, 1, 2, 4, 6, 5, 7 },
    /* 0b01011010 */ { 1, 3, 0, 2, 4, 6, 5, 7 },
    /* 0b01011011 */ { 0, 1, 6, 2, 3, 4, 5, 7 },
    /* 0b01011100 */ { 2, 3, 0, 1, 4, 6, 5, 7 },
    /* 0b01011101 */ { 0, 2, 6, 1, 3, 4, 5, 7 },
    /* 0b01011110 */ { 1, 2, 6, 0, 3, 4, 5, 7 },
    /* 0b01011111 */ { 0, 1, 4, 6, 2, 3, 5, 7 },
    /* 0b01100000 */ { 5, 6, 2, 3, 0, 1, 4, 7 },
    /* 0b01100001 */ { 0, 5, 2, 3, 6, 1, 4, 7 },
    /* 0b01100010 */ { 1, 5, 2, 3, 6, 0, 4, 7 },
    /* 0b01100011 */ { 0, 1, 2, 3, 5, 6, 4, 7 },
    /* 0b01100100 */ { 2, 5, 1, 3, 6, 0, 4, 7 },
    /* 0b01100101 */ { 0, 2, 1, 3, 5, 6, 4, 7 },
    /* 0b01100110 */ { 1, 2, 0, 3, 5, 6, 4, 7 },
    /* 0b01100111 */ { 0, 1, 6, 3, 2, 5, 4, 7 },
    /* 0b01101000 */ { 3, 5, 1, 2, 6, 0, 4, 7 },
    /* 0b01101001 */ { 0, 3, 1, 2, 5, 6, 4, 7 },
    /* 0b01101010 */ { 1, 3, 0, 2, 5, 6, 4, 7 },
    /* 0b01101011 */ { 0, 1, 6, 2, 3, 5, 4, 7 },
    /* 0b01101100 */ { 2, 3, 0, 1, 5, 6, 4, 7 },
    /* 0b01101101 */ { 0, 2, 6, 1, 3, 5, 4, 7 },
    /* 0b01101110 */ { 1, 2, 6, 0, 3, 5, 4, 7 },
    /* 0b01101111 */ { 0, 1, 5, 6, 2, 3, 4, 7 },
    /* 0b01110000 */ { 4, 5, 1, 2, 6, 0, 3, 7 },
    /* 0b01110001 */ { 0, 4, 1, 2, 5, 6, 3, 7 },
    /* 0b01110010 */ { 1, 4, 0, 2, 5, 6, 3, 7 },
    /* 0b01110011 */ { 0, 1, 6, 2, 4, 5, 3, 7 },
    /* 0b01110100 */ { 2, 4, 0, 1, 5, 6, 3, 7 },
    /* 0b01110101 */ { 0, 2, 6, 1, 4, 5, 3, 7 },
    /* 0b01110110 */ { 1, 2, 6, 0, 4, 5, 3, 7 },
    /* 0b01110111 */ { 0, 1, 5, 6, 2, 4, 3, 7 },
    /* 0b01111000 */ { 3, 4, 0, 1, 5, 6, 2, 7 },
    /* 0b01111001 */ { 0, 3, 6, 1, 4, 5, 2, 7 },
    /* 0b01111010 */ { 1, 3, 6, 0, 4, 5, 2, 7 },
    /* 0b01111011 */ { 0, 1, 5, 6, 3, 4, 2, 7 },
    /* 0b01111100 */ { 2, 3, 6, 0, 4, 5, 1, 7 },
    /* 0b01111101 */ { 0, 2, 5, 6, 3, 4, 1, 7 },
    /* 0b01111110 */ { 1, 2, 5, 6, 3, 4, 0, 7 },
    /* 0b01111111 */ { 0, 1, 4, 5, 2, 3, 6, 7 },
    /* 0b10000000 */ { 7, 0, 3, 4, 1, 2, 5, 6 },
    /* 0b10000001 */ { 0, 7, 3, 4, 1, 2, 5, 6 },
    /* 0b10000010 */ { 1, 7, 3, 4, 0, 2, 5, 6 },
    /* 0b10000011 */ { 0, 1, 3, 4, 7, 2, 5, 6 },
    /* 0b10000100 */ { 2, 7, 3, 4, 0, 1, 5, 6 },
    /* 0b10000101 */ { 0, 2, 3, 4, 7, 1, 5, 6 },
    /* 0b10000110 */ { 1, 2, 3, 4, 7, 0, 5, 6 },
    /* 0b10000111 */ { 0, 1, 3, 4, 2, 7, 5, 6 },
    /* 0b10001000 */ { 3, 7, 2, 4, 0, 1, 5, 6 },
    /* 0b10001001 */ { 0, 3, 2, 4, 7, 1, 5, 6 },
    /* 0b10001010 */ { 1, 3, 2, 4, 7, 0, 5, 6 },
    /* 0b10001011 */ { 0, 1, 2, 4, 3, 7, 5, 6 },
    /* 0b10001100 */ { 2, 3, 1, 4, 7, 0, 5, 6 },
    /* 0b10001101 */ { 0, 2, 1, 4, 3, 7, 5, 6 },
    /* 0b10001110 */ { 1, 2, 0, 4, 3, 7, 5, 6 },
    /* 0b10001111 */ { 0, 1, 7, 4, 2, 3, 5, 6 },
    /* 0b10010000 */ { 4, 7, 2, 3, 0, 1, 5, 6 },
    /* 0b10010001 */ { 0, 4, 2, 3, 7, 1, 5, 6 },
    /* 0b10010010 */ { 1, 4, 2, 3, 7, 0, 5, 6 },
    /* 0b10010011 */ { 0, 1, 2, 3, 4, 7, 5, 6 },
    /* 0b10010100 */ { 2, 4, 1, 3, 7, 0, 5, 6 },
    /* 0b10010101 */ { 0, 2, 1, 3, 4, 7, 5, 6 },
    /* 0b10010110 */ { 1, 2, 0, 3, 4, 7, 5, 6 },
    /* 0b10010111 */ { 0, 1, 7, 3, 2, 4, 5, 6 },
    /* 0b10011000 */ { 3, 4, 1, 2, 7, 0, 5, 6 },
    /* 0b10011001 */ { 0, 3, 1, 2, 4, 7, 5, 6 },
    /* 0b10011010 */ { 1, 3, 0, 2, 4, 7, 5, 6 },
    /* 0b10011011 */ { 0, 1, 7, 2, 3, 4, 5, 6 },
    /* 0b10011100 */ { 2, 3, 0, 1, 4, 7, 5, 6 },
    /* 0b10011101 */ { 0, 2, 7, 1, 3, 4, 5, 6 },
    /* 0b10011110 */ { 1, 2, 7, 0, 3, 4, 5, 6 },
    /* 0b10011111 */ { 0, 1, 4, 7, 2, 3, 5, 6 },
    /* 0b10100000 */ { 5, 7, 2, 3, 0, 1, 4, 6 },
    /* 0b10100001 */ { 0, 5, 2, 3, 7, 1, 4, 6 },
    /* 0b10100010 */ { 1, 5, 2, 3, 7, 0, 4, 6 },
    /* 0b10100011 */ { 0, 1, 2, 3, 5, 7, 4, 6 },
    /* 0b10100100 */ { 2, 5, 1, 3, 7, 0, 4, 6 },
    /* 0b10100101 */ { 0, 2, 1, 3, 5, 7, 4, 6 },
    /* 0b10100110 */ { 1, 2, 0, 3, 5, 7, 4, 6 },
    /* 0b10100111 */ { 0, 1, 7, 3, 2, 5, 4, 6 },
    /* 0b10101000 */ { 3, 5, 1, 2, 7, 0, 4, 6 },
    /* 0b10101001 */ { 0, 3, 1, 2, 5, 7, 4, 6 },
    /* 0b10101010 */ { 1, 3, 0, 2, 5, 7, 4, 6 },
    /* 0b10101011 */ { 0, 1, 7, 2, 3, 5, 4, 6 },
    /* 0b10101100 */ { 2, 3, 0, 1, 5, 7, 4, 6 },
    /* 0b10101101 */ { 0, 2, 7, 1, 3, 5, 4, 6 },
    /* 0b10101110 */ { 1, 2, 7, 0, 3, 5, 4, 6 },
    /* 0b10101111 */ { 0, 1, 5, 7, 2, 3, 4, 6 },
    /* 0b10110000 */ { 4, 5, 1, 2, 7, 0, 3, 6 },
    /* 0b10110001 */ { 0, 4, 1, 2, 5, 7, 3, 6 },
    /* 0b10110010 */ { 1, 4, 0, 2, 5, 7, 3, 6 },
    /* 0b10110011 */ { 0, 1, 7, 2, 4, 5, 3, 6 },
    /* 0b10110100 */ { 2, 4, 0, 1, 5, 7, 3, 6 },
    /* 0b10110101 */ { 0, 2, 7, 1, 4, 5, 3, 6 },
    /* 0b10110110 */ { 1, 2, 7, 0, 4, 5, 3, 6 },
    /* 0b10110111 */ { 0, 1, 5, 7, 2, 4, 3, 6 },
    /* 0b10111000 */ { 3, 4, 0, 1, 5, 7, 2, 6 },
    /* 0b10111001 */ { 0, 3, 7, 1, 4, 5, 2, 6 },
    /* 0b10111010 */ { 1, 3, 7, 0, 4, 5, 2, 6 },
    /* 0b10111011 */ { 0, 1, 5, 7, 3, 4, 2, 6 },
    /* 0b10111100 */ { 2, 3, 7, 0, 4, 5, 1, 6 },
    /* 0b10111101 */ { 0, 2, 5, 7, 3, 4, 1, 6 },
    /* 0b10111110 */ { 1, 2, 5, 7, 3, 4, 0, 6 },
    /* 0b10111111 */ { 0, 1, 4, 5, 2, 3, 7, 6 },
    /* 0b11000000 */ { 6, 7, 2, 3, 0, 1, 4, 5 },
    /* 0b11000001 */ { 0, 6, 2, 3, 7, 1, 4, 5 },
    /* 0b11000010 */ { 1, 6, 2, 3, 7, 0, 4, 5 },
    /* 0b11000011 */ { 0, 1, 2, 3, 6, 7, 4, 5 },
    /* 0b11000100 */ { 2, 6, 1, 3, 7, 0, 4, 5 },
    /* 0b11000101 */ { 0, 2, 1, 3, 6, 7, 4, 5 },
    /* 0b11000110 */ { 1, 2, 0, 3, 6, 7, 4, 5 },
    /* 0b11000111 */ { 0, 1, 7, 3, 2, 6, 4, 5 },
    /* 0b11001000 */ { 3, 6, 1, 2, 7, 0, 4, 5 },
    /* 0b11001001 */ { 0, 3, 1, 2, 6, 7, 4, 5 },
    /* 0b11001010 */ { 1, 3, 0, 2, 6, 7, 4, 5 },
    /* 0b11001011 */ { 0, 1, 7, 2, 3, 6, 4, 5 },
    /* 0b11001100 */ { 2, 3, 0, 1, 6, 7, 4, 5 },
    /* 0b11001101 */ { 0, 2, 7, 1, 3, 6, 4, 5 },
    /* 0b11001110 */ { 1, 2, 7, 0, 3, 6, 4, 5 },
    /* 0b11001111 */ { 0, 1, 6, 7, 2, 3, 4, 5 },
    /* 0b11010000 */ { 4, 6, 1, 2, 7, 0, 3, 5 },
    /* 0b11010001 */ { 0, 4, 1, 2, 6, 7, 3, 5 },
    /* 0b11010010 */ { 1, 4, 0, 2, 6, 7, 3, 5 },
    /* 0b11010011 */ { 0, 1, 7, 2, 4, 6, 3, 5 },
    /* 0b11010100 */ { 2, 4, 0, 1, 6, 7, 3, 5 },
    /* 0b11010101 */ { 0, 2, 7, 1, 4, 6, 3, 5 },
    /* 0b11010110 */ { 1, 2, 7, 0, 4, 6, 3, 5 },
    /* 0b11010111 */ { 0, 1, 6, 7, 2, 4, 3, 5 },
    /* 0b11011000 */ { 3, 4, 0, 1, 6, 7, 2, 5 },
    /* 0b11011001 */ { 0, 3, 7, 1, 4, 6, 2, 5 },
    /* 0b11011010 */ { 1, 3, 7, 0, 4, 6, 2, 5 },
    /* 0b11011011 */ { 0, 1, 6, 7, 3, 4, 2, 5 },
    /* 0b11011100 */ { 2, 3, 7, 0, 4, 6, 1, 5 },
    /* 0b11011101 */ { 0, 2, 6, 7, 3, 4, 1, 5 },
    /* 0b11011110 */ { 1, 2, 6, 7, 3, 4, 0, 5 },
    /* 0b11011111 */ { 0, 1, 4, 6, 2, 3, 7, 5 },
    /* 0b11100000 */ { 5, 6, 1, 2, 7, 0, 3, 4 },
    /* 0b11100001 */ { 0, 5, 1, 2, 6, 7, 3, 4 },
    /* 0b11100010 */ { 1, 5, 0, 2, 6, 7, 3, 4 },
    /* 0b11100011 */ { 0, 1, 7, 2, 5, 6, 3, 4 },
    /* 0b11100100 */ { 2, 5, 0, 1, 6, 7, 3, 4 },
    /* 0b11100101 */ { 0, 2, 7, 1, 5, 6, 3, 4 },
    /* 0b11100110 */ { 1, 2, 7, 0, 5, 6, 3, 4 },
    /* 0b11100111 */ { 0, 1, 6, 7, 2, 5, 3, 4 },
    /* 0b11101000 */ { 3, 5, 0, 1, 6, 7, 2, 4 },
    /* 0b11101001 */ { 0, 3, 7, 1, 5, 6, 2, 4 },
    /* 0b11101010 */ { 1, 3, 7, 0, 5, 6, 2, 4 },
    /* 0b11101011 */ { 0, 1, 6, 7, 3, 5, 2, 4 },
    /* 0b11101100 */ { 2, 3, 7, 0, 5, 6, 1, 4 },
    /* 0b11101101 */ { 0, 2, 6, 7, 3, 5, 1, 4 },
    /* 0b11101110 */ { 1, 2, 6, 7, 3, 5, 0, 4 },
    /* 0b11101111 */ { 0, 1, 5, 6, 2, 3, 7, 4 },
    /* 0b11110000 */ { 4, 5, 0, 1, 6, 7, 2, 3 },
    /* 0b11110001 */ { 0, 4, 7, 1, 5, 6, 2, 3 },
    /* 0b11110010 */ { 1, 4, 7, 0, 5, 6, 2, 3 },
    /* 0b11110011 */ { 0, 1, 6, 7, 4, 5, 2, 3 },
    /* 0b11110100 */ { 2, 4, 7, 0, 5, 6, 1, 3 },
    /* 0b11110101 */ { 0, 2, 6, 7, 4, 5, 1, 3 },
    /* 0b11110110 */ { 1, 2, 6, 7, 4, 5, 0, 3 },
    /* 0b11110111 */ { 0, 1, 5, 6, 2, 4, 7, 3 },
    /* 0b11111000 */ { 3, 4, 7, 0, 5, 6, 1, 2 },
    /* 0b11111001 */ { 0, 3, 6, 7, 4, 5, 1, 2 },
    /* 0b11111010 */ { 1, 3, 6, 7, 4, 5, 0, 2 },
    /* 0b11111011 */ { 0, 1, 5, 6, 3, 4, 7, 2 },
    /* 0b11111100 */ { 2, 3, 6, 7, 4, 5, 0, 1 },
    /* 0b11111101 */ { 0, 2, 5, 6, 3, 4, 7, 1 },
    /* 0b11111110 */ { 1, 2, 5, 6, 3, 4, 7, 0 },
    /* 0b11111111 */ { 0, 1, 4, 5, 2, 3, 6, 7 },
};
