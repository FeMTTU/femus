/* /////////////////////////////////////////////////////////////////////////
 * File:        test.unit.cpp.1.cpp
 *
 * Purpose:     Implementation file for the test.unit.cpp.1 project.
 *
 * Created:     21st April 2009
 * Updated:     20th January 2010
 *
 * Status:      Wizard-generated
 *
 * License:     (Licensed under the Synesis Software Open License)
 *
 *              Copyright (c) 2009-2010, Synesis Software Pty Ltd.
 *              All rights reserved.
 *
 *              www:        http://www.synesis.com.au/software
 *
 * ////////////////////////////////////////////////////////////////////// */


/* /////////////////////////////////////////////////////////////////////////
 * Test component header file include(s)
 */

#include <b64/b64.hpp>

/* /////////////////////////////////////////////////////////////////////////
 * Includes
 */

/* xTests Header Files */
#include <xtests/xtests.h>

/* STLSoft Header Files */
#include <stlsoft/stlsoft.h>
#include <stlsoft/string/simple_string.hpp>

/* Standard C Header Files */
#include <stdlib.h>

/* /////////////////////////////////////////////////////////////////////////
 * Forward declarations
 */

namespace 
{

static void test_encode_NULL();
static void test_encode_empty();
static void test_encode_1_0();
static void test_encode_1_1();
static void test_encode_5_1();
static void test_encode_5_2();
static void test_encode_8_1();
static void test_encode_8_2();
static void test_encode_8_3();
static void test_encode_8_4();
static void test_encode_80_1();
static void test_encode_80_2();
static void test_encode_80_3();
static void test_encode_80_4();

static void test_1_0(void);
static void test_1_1(void);
static void test_1_2(void);
static void test_1_3(void);
static void test_1_4(void);
static void test_1_5(void);
static void test_1_6(void);
static void test_1_7(void);
static void test_1_8(void);
static void test_1_9(void);

static void test_2_0(void);
static void test_2_1(void);
static void test_2_2(void);
static void test_2_3(void);
static void test_2_4(void);
static void test_2_5(void);
static void test_2_6(void);
static void test_2_7(void);
static void test_2_8(void);
static void test_2_9(void);

static void test_3_00(void);
static void test_3_01(void);
static void test_3_02(void);
static void test_3_03(void);
static void test_3_04(void);
static void test_3_05(void);
static void test_3_06(void);
static void test_3_07(void);
static void test_3_08(void);
static void test_3_09(void);
static void test_3_10(void);

#if 0
static void test_4_0(void);
static void test_4_1(void);
static void test_4_2(void);
static void test_4_3(void);
static void test_4_4(void);
static void test_4_5(void);
static void test_4_6(void);
static void test_4_7(void);
static void test_4_8(void);
static void test_4_9(void);
#endif /* 0 */

} /* anonymous namespace */

/* /////////////////////////////////////////////////////////////////////////
 * Main
 */

int main(int argc, char **argv)
{
    int retCode = EXIT_SUCCESS;
    int verbosity = 2;

    XTESTS_COMMANDLINE_PARSEVERBOSITY(argc, argv, &verbosity);

    if(XTESTS_START_RUNNER("test.unit.cpp.1", verbosity))
    {
        XTESTS_RUN_CASE(test_encode_NULL);
        XTESTS_RUN_CASE(test_encode_empty);
        XTESTS_RUN_CASE(test_encode_1_0);
        XTESTS_RUN_CASE(test_encode_1_1);
        XTESTS_RUN_CASE(test_encode_5_1);
        XTESTS_RUN_CASE(test_encode_5_2);
        XTESTS_RUN_CASE(test_encode_8_1);
        XTESTS_RUN_CASE(test_encode_8_2);
        XTESTS_RUN_CASE(test_encode_8_3);
        XTESTS_RUN_CASE(test_encode_8_4);
        XTESTS_RUN_CASE(test_encode_80_1);
        XTESTS_RUN_CASE(test_encode_80_2);
        XTESTS_RUN_CASE(test_encode_80_3);
        XTESTS_RUN_CASE(test_encode_80_4);


        XTESTS_RUN_CASE(test_1_0);
        XTESTS_RUN_CASE(test_1_1);
        XTESTS_RUN_CASE(test_1_2);
        XTESTS_RUN_CASE(test_1_3);
        XTESTS_RUN_CASE(test_1_4);
        XTESTS_RUN_CASE(test_1_5);
        XTESTS_RUN_CASE(test_1_6);
        XTESTS_RUN_CASE(test_1_7);
        XTESTS_RUN_CASE(test_1_8);
        XTESTS_RUN_CASE(test_1_9);

        XTESTS_RUN_CASE(test_2_0);
        XTESTS_RUN_CASE(test_2_1);
        XTESTS_RUN_CASE(test_2_2);
        XTESTS_RUN_CASE(test_2_3);
        XTESTS_RUN_CASE(test_2_4);
        XTESTS_RUN_CASE(test_2_5);
        XTESTS_RUN_CASE(test_2_6);
        XTESTS_RUN_CASE(test_2_7);
        XTESTS_RUN_CASE(test_2_8);
        XTESTS_RUN_CASE(test_2_9);

        XTESTS_RUN_CASE(test_3_00);
        XTESTS_RUN_CASE(test_3_01);
        XTESTS_RUN_CASE(test_3_02);
        XTESTS_RUN_CASE(test_3_03);
        XTESTS_RUN_CASE(test_3_04);
        XTESTS_RUN_CASE(test_3_05);
        XTESTS_RUN_CASE(test_3_06);
        XTESTS_RUN_CASE(test_3_07);
        XTESTS_RUN_CASE(test_3_08);
        XTESTS_RUN_CASE(test_3_09);
        XTESTS_RUN_CASE(test_3_10);

#if 0
        XTESTS_RUN_CASE(test_4_0);
        XTESTS_RUN_CASE(test_4_1);
        XTESTS_RUN_CASE(test_4_2);
        XTESTS_RUN_CASE(test_4_3);
        XTESTS_RUN_CASE(test_4_4);
        XTESTS_RUN_CASE(test_4_5);
        XTESTS_RUN_CASE(test_4_6);
        XTESTS_RUN_CASE(test_4_7);
        XTESTS_RUN_CASE(test_4_8);
        XTESTS_RUN_CASE(test_4_9);
#endif /* 0 */

        XTESTS_PRINT_RESULTS();

        XTESTS_END_RUNNER_UPDATE_EXITCODE(&retCode);
    }

    return retCode;
}

/* /////////////////////////////////////////////////////////////////////////
 * Test function implementations
 */

namespace 
{

//  inline string_t encode(void const *src, size_t srcSize);
//
//  inline string_t encode(void const *src, size_t srcSize, unsigned flags, int lineLen = 0, B64_RC *rc = NULL);
//
//  template <typename T, size_t N>
//  inline string_t encode(T (&ar)[N]);
//
//  inline string_t encode(blob_t const &blob);
//
//  inline string_t encode(blob_t const &blob, unsigned flags, int lineLen = 0, B64_RC *rc = NULL);


static void test_encode_NULL()
{
    XTESTS_TEST_MULTIBYTE_STRING_EQUAL("", b64::encode(NULL, 0u));

    XTESTS_TEST_MULTIBYTE_STRING_EQUAL("", b64::encode(NULL, 0u, 0));
    XTESTS_TEST_MULTIBYTE_STRING_EQUAL("", b64::encode(NULL, 0u, 0, -1));
    XTESTS_TEST_MULTIBYTE_STRING_EQUAL("", b64::encode(NULL, 0u, 0, -1, NULL));
}

static void test_encode_empty()
{
    unsigned int    bytes[1];

    XTESTS_TEST_MULTIBYTE_STRING_EQUAL("", b64::encode(&bytes[0], 0u));

    XTESTS_TEST_MULTIBYTE_STRING_EQUAL("", b64::encode(&bytes[0], 0u, 0));
    XTESTS_TEST_MULTIBYTE_STRING_EQUAL("", b64::encode(&bytes[0], 0u, 0, -1));
    XTESTS_TEST_MULTIBYTE_STRING_EQUAL("", b64::encode(&bytes[0], 0u, 0, -1, NULL));

    b64::blob_t blob;

    XTESTS_TEST_MULTIBYTE_STRING_EQUAL("", b64::encode(blob));

    XTESTS_TEST_MULTIBYTE_STRING_EQUAL("", b64::encode(blob, 0));
    XTESTS_TEST_MULTIBYTE_STRING_EQUAL("", b64::encode(blob, 0, -1));
    XTESTS_TEST_MULTIBYTE_STRING_EQUAL("", b64::encode(blob, 0, -1, NULL));
}

static void test_encode_1_0()
{
    unsigned char   bytes[] = { 0 };

    XTESTS_TEST_MULTIBYTE_STRING_EQUAL("AA==", b64::encode(&bytes[0], sizeof(bytes)));

    XTESTS_TEST_MULTIBYTE_STRING_EQUAL("AA==", b64::encode(&bytes[0], sizeof(bytes), 0));
    XTESTS_TEST_MULTIBYTE_STRING_EQUAL("AA==", b64::encode(&bytes[0], sizeof(bytes), 0, -1));
    XTESTS_TEST_MULTIBYTE_STRING_EQUAL("AA==", b64::encode(&bytes[0], sizeof(bytes), 0, -1, NULL));

    b64::blob_t blob(STLSOFT_NUM_ELEMENTS(bytes));

    std::copy(&bytes[0], &bytes[0] + STLSOFT_NUM_ELEMENTS(bytes), blob.begin());

    XTESTS_TEST_MULTIBYTE_STRING_EQUAL("AA==", b64::encode(blob));

    XTESTS_TEST_MULTIBYTE_STRING_EQUAL("AA==", b64::encode(blob, 0));
    XTESTS_TEST_MULTIBYTE_STRING_EQUAL("AA==", b64::encode(blob, 0, -1));
    XTESTS_TEST_MULTIBYTE_STRING_EQUAL("AA==", b64::encode(blob, 0, -1, NULL));
}

static void test_encode_1_1()
{
    unsigned char   bytes[] = { 1 };

    XTESTS_TEST_MULTIBYTE_STRING_EQUAL("AQ==", b64::encode(&bytes[0], sizeof(bytes)));

    XTESTS_TEST_MULTIBYTE_STRING_EQUAL("AQ==", b64::encode(&bytes[0], sizeof(bytes), 0));
    XTESTS_TEST_MULTIBYTE_STRING_EQUAL("AQ==", b64::encode(&bytes[0], sizeof(bytes), 0, -1));
    XTESTS_TEST_MULTIBYTE_STRING_EQUAL("AQ==", b64::encode(&bytes[0], sizeof(bytes), 0, -1, NULL));

    b64::blob_t blob(STLSOFT_NUM_ELEMENTS(bytes));

    std::copy(&bytes[0], &bytes[0] + STLSOFT_NUM_ELEMENTS(bytes), blob.begin());

    XTESTS_TEST_MULTIBYTE_STRING_EQUAL("AQ==", b64::encode(blob));

    XTESTS_TEST_MULTIBYTE_STRING_EQUAL("AQ==", b64::encode(blob, 0));
    XTESTS_TEST_MULTIBYTE_STRING_EQUAL("AQ==", b64::encode(blob, 0, -1));
    XTESTS_TEST_MULTIBYTE_STRING_EQUAL("AQ==", b64::encode(blob, 0, -1, NULL));
}

static void test_encode_5_1()
{
    unsigned char   bytes[] = { 1, 1, 1, 1, 1 };

    XTESTS_TEST_MULTIBYTE_STRING_EQUAL("AQEBAQE=", b64::encode(&bytes[0], sizeof(bytes)));

    XTESTS_TEST_MULTIBYTE_STRING_EQUAL("AQEBAQE=", b64::encode(&bytes[0], sizeof(bytes), 0));
    XTESTS_TEST_MULTIBYTE_STRING_EQUAL("AQEBAQE=", b64::encode(&bytes[0], sizeof(bytes), 0, -1));
    XTESTS_TEST_MULTIBYTE_STRING_EQUAL("AQEBAQE=", b64::encode(&bytes[0], sizeof(bytes), 0, -1, NULL));

    b64::blob_t blob(STLSOFT_NUM_ELEMENTS(bytes));

    std::copy(&bytes[0], &bytes[0] + STLSOFT_NUM_ELEMENTS(bytes), blob.begin());

    XTESTS_TEST_MULTIBYTE_STRING_EQUAL("AQEBAQE=", b64::encode(blob));

    XTESTS_TEST_MULTIBYTE_STRING_EQUAL("AQEBAQE=", b64::encode(blob, 0));
    XTESTS_TEST_MULTIBYTE_STRING_EQUAL("AQEBAQE=", b64::encode(blob, 0, -1));
    XTESTS_TEST_MULTIBYTE_STRING_EQUAL("AQEBAQE=", b64::encode(blob, 0, -1, NULL));
}

static void test_encode_5_2()
{
    unsigned char   bytes[] = { 1, 2, 2, 2, 2 };

    XTESTS_TEST_MULTIBYTE_STRING_EQUAL("AQICAgI=", b64::encode(&bytes[0], sizeof(bytes)));

    XTESTS_TEST_MULTIBYTE_STRING_EQUAL("AQICAgI=", b64::encode(&bytes[0], sizeof(bytes), 0));
    XTESTS_TEST_MULTIBYTE_STRING_EQUAL("AQICAgI=", b64::encode(&bytes[0], sizeof(bytes), 0, -1));
    XTESTS_TEST_MULTIBYTE_STRING_EQUAL("AQICAgI=", b64::encode(&bytes[0], sizeof(bytes), 0, -1, NULL));

    b64::blob_t blob(STLSOFT_NUM_ELEMENTS(bytes));

    std::copy(&bytes[0], &bytes[0] + STLSOFT_NUM_ELEMENTS(bytes), blob.begin());

    XTESTS_TEST_MULTIBYTE_STRING_EQUAL("AQICAgI=", b64::encode(blob));

    XTESTS_TEST_MULTIBYTE_STRING_EQUAL("AQICAgI=", b64::encode(blob, 0));
    XTESTS_TEST_MULTIBYTE_STRING_EQUAL("AQICAgI=", b64::encode(blob, 0, -1));
    XTESTS_TEST_MULTIBYTE_STRING_EQUAL("AQICAgI=", b64::encode(blob, 0, -1, NULL));
}

static void test_encode_8_1()
{
    unsigned char   bytes[] = { 1, 1, 1, 1, 1, 1, 1, 1 };

    XTESTS_TEST_MULTIBYTE_STRING_EQUAL("AQEBAQEBAQE=", b64::encode(&bytes[0], sizeof(bytes)));

    XTESTS_TEST_MULTIBYTE_STRING_EQUAL("AQEBAQEBAQE=", b64::encode(&bytes[0], sizeof(bytes), 0));
    XTESTS_TEST_MULTIBYTE_STRING_EQUAL("AQEBAQEBAQE=", b64::encode(&bytes[0], sizeof(bytes), 0, -1));
    XTESTS_TEST_MULTIBYTE_STRING_EQUAL("AQEBAQEBAQE=", b64::encode(&bytes[0], sizeof(bytes), 0, -1, NULL));

    b64::blob_t blob(STLSOFT_NUM_ELEMENTS(bytes));

    std::copy(&bytes[0], &bytes[0] + STLSOFT_NUM_ELEMENTS(bytes), blob.begin());

    XTESTS_TEST_MULTIBYTE_STRING_EQUAL("AQEBAQEBAQE=", b64::encode(blob));

    XTESTS_TEST_MULTIBYTE_STRING_EQUAL("AQEBAQEBAQE=", b64::encode(blob, 0));
    XTESTS_TEST_MULTIBYTE_STRING_EQUAL("AQEBAQEBAQE=", b64::encode(blob, 0, -1));
    XTESTS_TEST_MULTIBYTE_STRING_EQUAL("AQEBAQEBAQE=", b64::encode(blob, 0, -1, NULL));
}

static void test_encode_8_2()
{
    unsigned char   bytes[] = { 1, 2, 3, 4, 5, 6, 7, 8 };

    XTESTS_TEST_MULTIBYTE_STRING_EQUAL("AQIDBAUGBwg=", b64::encode(&bytes[0], sizeof(bytes)));

    XTESTS_TEST_MULTIBYTE_STRING_EQUAL("AQIDBAUGBwg=", b64::encode(&bytes[0], sizeof(bytes), 0));
    XTESTS_TEST_MULTIBYTE_STRING_EQUAL("AQIDBAUGBwg=", b64::encode(&bytes[0], sizeof(bytes), 0, -1));
    XTESTS_TEST_MULTIBYTE_STRING_EQUAL("AQIDBAUGBwg=", b64::encode(&bytes[0], sizeof(bytes), 0, -1, NULL));

    b64::blob_t blob(STLSOFT_NUM_ELEMENTS(bytes));

    std::copy(&bytes[0], &bytes[0] + STLSOFT_NUM_ELEMENTS(bytes), blob.begin());

    XTESTS_TEST_MULTIBYTE_STRING_EQUAL("AQIDBAUGBwg=", b64::encode(blob));

    XTESTS_TEST_MULTIBYTE_STRING_EQUAL("AQIDBAUGBwg=", b64::encode(blob, 0));
    XTESTS_TEST_MULTIBYTE_STRING_EQUAL("AQIDBAUGBwg=", b64::encode(blob, 0, -1));
    XTESTS_TEST_MULTIBYTE_STRING_EQUAL("AQIDBAUGBwg=", b64::encode(blob, 0, -1, NULL));
}

static void test_encode_8_3()
{
    unsigned char   bytes[] = { 1, 2, 3, 4, 5, 6, 7, 8 };

    XTESTS_TEST_MULTIBYTE_STRING_EQUAL("AQIDBAUGBwg=", b64::encode(&bytes[0], sizeof(bytes), 0, 40));
    XTESTS_TEST_MULTIBYTE_STRING_EQUAL("AQIDBAUGBwg=", b64::encode(&bytes[0], sizeof(bytes), 0, 40, NULL));

    b64::blob_t blob(STLSOFT_NUM_ELEMENTS(bytes));

    std::copy(&bytes[0], &bytes[0] + STLSOFT_NUM_ELEMENTS(bytes), blob.begin());

    XTESTS_TEST_MULTIBYTE_STRING_EQUAL("AQIDBAUGBwg=", b64::encode(blob, 0, 40));
    XTESTS_TEST_MULTIBYTE_STRING_EQUAL("AQIDBAUGBwg=", b64::encode(blob, 0, 40, NULL));
}

static void test_encode_8_4()
{
    unsigned char   bytes[] = { 1, 2, 3, 4, 5, 6, 7, 8 };

    XTESTS_TEST_MULTIBYTE_STRING_EQUAL("AQID\r\nBAUG\r\nBwg=", b64::encode(&bytes[0], sizeof(bytes), 0, 4));
    XTESTS_TEST_MULTIBYTE_STRING_EQUAL("AQID\r\nBAUG\r\nBwg=", b64::encode(&bytes[0], sizeof(bytes), 0, 4, NULL));

    b64::blob_t blob(STLSOFT_NUM_ELEMENTS(bytes));

    std::copy(&bytes[0], &bytes[0] + STLSOFT_NUM_ELEMENTS(bytes), blob.begin());

    XTESTS_TEST_MULTIBYTE_STRING_EQUAL("AQID\r\nBAUG\r\nBwg=", b64::encode(blob, 0, 4));
    XTESTS_TEST_MULTIBYTE_STRING_EQUAL("AQID\r\nBAUG\r\nBwg=", b64::encode(blob, 0, 4, NULL));
}

static void test_encode_80_1()
{
    unsigned char   bytes[] = { 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1 };

    XTESTS_TEST_MULTIBYTE_STRING_EQUAL("AQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQE=", b64::encode(&bytes[0], sizeof(bytes)));

    XTESTS_TEST_MULTIBYTE_STRING_EQUAL("AQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQE=", b64::encode(&bytes[0], sizeof(bytes), 0));
    XTESTS_TEST_MULTIBYTE_STRING_EQUAL("AQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEB\r\nAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQE=", b64::encode(&bytes[0], sizeof(bytes), 0, -1));
    XTESTS_TEST_MULTIBYTE_STRING_EQUAL("AQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEB\r\nAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQE=", b64::encode(&bytes[0], sizeof(bytes), 0, -1, NULL));

    b64::blob_t blob(STLSOFT_NUM_ELEMENTS(bytes));

    std::copy(&bytes[0], &bytes[0] + STLSOFT_NUM_ELEMENTS(bytes), blob.begin());

    XTESTS_TEST_MULTIBYTE_STRING_EQUAL("AQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQE=", b64::encode(blob));

    XTESTS_TEST_MULTIBYTE_STRING_EQUAL("AQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQE=", b64::encode(blob, 0));
    XTESTS_TEST_MULTIBYTE_STRING_EQUAL("AQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEB\r\nAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQE=", b64::encode(blob, 0, -1));
    XTESTS_TEST_MULTIBYTE_STRING_EQUAL("AQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEB\r\nAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQE=", b64::encode(blob, 0, -1, NULL));
}

static void test_encode_80_2()
{
    unsigned char   bytes[] = { 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1 };

    XTESTS_TEST_MULTIBYTE_STRING_EQUAL("AQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQE=", b64::encode(&bytes[0], sizeof(bytes)));

    XTESTS_TEST_MULTIBYTE_STRING_EQUAL("AQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEB\r\nAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQE=", b64::encode(&bytes[0], sizeof(bytes), b64::B64_F_LINE_LEN_64));
    XTESTS_TEST_MULTIBYTE_STRING_EQUAL("AQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEB\r\nAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQE=", b64::encode(&bytes[0], sizeof(bytes), b64::B64_F_LINE_LEN_64, 0));
    XTESTS_TEST_MULTIBYTE_STRING_EQUAL("AQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEB\r\nAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQE=", b64::encode(&bytes[0], sizeof(bytes), b64::B64_F_LINE_LEN_64, 0, NULL));
    XTESTS_TEST_MULTIBYTE_STRING_EQUAL("AQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEB\r\nAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQE=", b64::encode(&bytes[0], sizeof(bytes), 0, 64));
    XTESTS_TEST_MULTIBYTE_STRING_EQUAL("AQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEB\r\nAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQE=", b64::encode(&bytes[0], sizeof(bytes), 0, 64, NULL));

    b64::blob_t blob(STLSOFT_NUM_ELEMENTS(bytes));

    std::copy(&bytes[0], &bytes[0] + STLSOFT_NUM_ELEMENTS(bytes), blob.begin());

    XTESTS_TEST_MULTIBYTE_STRING_EQUAL("AQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQE=", b64::encode(blob));

    XTESTS_TEST_MULTIBYTE_STRING_EQUAL("AQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEB\r\nAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQE=", b64::encode(blob, b64::B64_F_LINE_LEN_64));
    XTESTS_TEST_MULTIBYTE_STRING_EQUAL("AQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEB\r\nAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQE=", b64::encode(blob, b64::B64_F_LINE_LEN_64, 0));
    XTESTS_TEST_MULTIBYTE_STRING_EQUAL("AQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEB\r\nAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQE=", b64::encode(blob, b64::B64_F_LINE_LEN_64, 0, NULL));
    XTESTS_TEST_MULTIBYTE_STRING_EQUAL("AQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEB\r\nAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQE=", b64::encode(blob, 0, 64));
    XTESTS_TEST_MULTIBYTE_STRING_EQUAL("AQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEB\r\nAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQE=", b64::encode(blob, 0, 64, NULL));
}

static void test_encode_80_3()
{
    unsigned char   bytes[] = { 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1 };

    XTESTS_TEST_MULTIBYTE_STRING_EQUAL("AQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQE=", b64::encode(&bytes[0], sizeof(bytes)));

    XTESTS_TEST_MULTIBYTE_STRING_EQUAL("AQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEB\r\nAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQE=", b64::encode(&bytes[0], sizeof(bytes), b64::B64_F_LINE_LEN_76));
    XTESTS_TEST_MULTIBYTE_STRING_EQUAL("AQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEB\r\nAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQE=", b64::encode(&bytes[0], sizeof(bytes), b64::B64_F_LINE_LEN_76, 0));
    XTESTS_TEST_MULTIBYTE_STRING_EQUAL("AQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEB\r\nAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQE=", b64::encode(&bytes[0], sizeof(bytes), b64::B64_F_LINE_LEN_76, 0, NULL));
    XTESTS_TEST_MULTIBYTE_STRING_EQUAL("AQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEB\r\nAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQE=", b64::encode(&bytes[0], sizeof(bytes), 0, 76));
    XTESTS_TEST_MULTIBYTE_STRING_EQUAL("AQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEB\r\nAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQE=", b64::encode(&bytes[0], sizeof(bytes), 0, 76, NULL));

    b64::blob_t blob(STLSOFT_NUM_ELEMENTS(bytes));

    std::copy(&bytes[0], &bytes[0] + STLSOFT_NUM_ELEMENTS(bytes), blob.begin());

    XTESTS_TEST_MULTIBYTE_STRING_EQUAL("AQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQE=", b64::encode(blob));

    XTESTS_TEST_MULTIBYTE_STRING_EQUAL("AQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEB\r\nAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQE=", b64::encode(blob, b64::B64_F_LINE_LEN_76));
    XTESTS_TEST_MULTIBYTE_STRING_EQUAL("AQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEB\r\nAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQE=", b64::encode(blob, 0, 76));
    XTESTS_TEST_MULTIBYTE_STRING_EQUAL("AQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEB\r\nAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQE=", b64::encode(blob, 0, 76, NULL));
}

static void test_encode_80_4()
{
    unsigned char   bytes[] = { 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1 };

    XTESTS_TEST_MULTIBYTE_STRING_EQUAL("AQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEB\r\nAQEBAQEBAQEBAQEBAQEBAQEBAQE=", b64::encode(&bytes[0], sizeof(bytes), 0, 80));
    XTESTS_TEST_MULTIBYTE_STRING_EQUAL("AQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEB\r\nAQEBAQEBAQEBAQEBAQEBAQEBAQE=", b64::encode(&bytes[0], sizeof(bytes), 0, 80, NULL));

    b64::blob_t blob(STLSOFT_NUM_ELEMENTS(bytes));

    std::copy(&bytes[0], &bytes[0] + STLSOFT_NUM_ELEMENTS(bytes), blob.begin());

    XTESTS_TEST_MULTIBYTE_STRING_EQUAL("AQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEB\r\nAQEBAQEBAQEBAQEBAQEBAQEBAQE=", b64::encode(blob, 0, 80));
    XTESTS_TEST_MULTIBYTE_STRING_EQUAL("AQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEB\r\nAQEBAQEBAQEBAQEBAQEBAQEBAQE=", b64::encode(blob, 0, 80, NULL));
}




static void test_1_0()
{
    XTESTS_TEST_MULTIBYTE_STRING_EQUAL("", b64::encode(NULL, 0u));
}

static void test_1_1()
{

//  XTESTS_TEST_MULTIBYTE_STRING_EQUAL("", b64::encode(&bytes[0], sizeof(bytes)));
}

static void test_1_2()
{
    unsigned char   bytes[] = { 1 };

    XTESTS_TEST_MULTIBYTE_STRING_EQUAL("AQ==", b64::encode(&bytes[0], sizeof(bytes)));
}

static void test_1_3()
{
    unsigned char   bytes[] = { 10 };

    XTESTS_TEST_MULTIBYTE_STRING_EQUAL("Cg==", b64::encode(&bytes[0], sizeof(bytes)));
}

static void test_1_4()
{
    unsigned char   bytes[] = { 100 };

    XTESTS_TEST_MULTIBYTE_STRING_EQUAL("ZA==", b64::encode(&bytes[0], sizeof(bytes)));
}

static void test_1_5()
{
    unsigned char   bytes[] = { 128 };

    XTESTS_TEST_MULTIBYTE_STRING_EQUAL("gA==", b64::encode(&bytes[0], sizeof(bytes)));
}

static void test_1_6()
{
    unsigned char   bytes[] = { 255 };

    XTESTS_TEST_MULTIBYTE_STRING_EQUAL("/w==", b64::encode(&bytes[0], sizeof(bytes)));
}

static void test_1_7()
{
    char bytes[] = "This is a test string";

    XTESTS_TEST_MULTIBYTE_STRING_EQUAL("VGhpcyBpcyBhIHRlc3Qgc3RyaW5n", b64::encode(&bytes[0], sizeof(bytes) - 1));
}

static void test_1_8()
{
    signed char bytes[] = { 7, 6, 5, 4, 3, 2, 1, 0, -1, -2, -3, -4, -5, -6, -7 };

    XTESTS_TEST_MULTIBYTE_STRING_EQUAL("BwYFBAMCAQD//v38+/r5", b64::encode(&bytes[0], sizeof(bytes)));
}

static void test_1_9()
{
}




static void test_2_0()
{
}

static void test_2_1()
{
}

static void test_2_2()
{
}

static void test_2_3()
{
}

static void test_2_4()
{
}

static void test_2_5()
{
}

static void test_2_6()
{
}

static void test_2_7()
{
}

static void test_2_8()
{
}

static void test_2_9()
{
}



static void test_3_00(void)
{
    using namespace b64;

    char const enc[] = "abcdefg%";

    try
    {
        blob_t r = decode(enc, STLSOFT_NUM_ELEMENTS(enc) - 1u);

        XTESTS_TEST_FAIL("should not reach here");
    }
    catch(b64::coding_exception& x)
    {
        XTESTS_TEST_ENUM_EQUAL(b64::B64_RC_DATA_ERROR, x.get_rc());
        XTESTS_TEST_POINTER_NOT_EQUAL(NULL, x.get_badChar());
        XTESTS_TEST_CHARACTER_EQUAL('%', *x.get_badChar());
    }
}

static void test_3_01(void)
{
    using namespace b64;

    char const enc[] = "abcdefg%";

    try
    {
        blob_t r = decode(enc);

        XTESTS_TEST_FAIL("should not reach here");
    }
    catch(b64::coding_exception& x)
    {
        XTESTS_TEST_ENUM_EQUAL(b64::B64_RC_DATA_ERROR, x.get_rc());
        XTESTS_TEST_POINTER_EQUAL(NULL, x.get_badChar());
    }
}

static void test_3_02(void)
{
    using namespace b64;

    char const enc[] = "abcdefg%";

    try
    {
        blob_t r = decode(B64_F_STOP_ON_NOTHING, enc);

        XTESTS_TEST_PASSED();
    }
    catch(b64::coding_exception& /* x */)
    {
        XTESTS_TEST_FAIL("should not reach here");
    }
}

static void test_3_03(void)
{
    using namespace b64;

    stlsoft::simple_string const enc = "abcdefg%";

    try
    {
        blob_t r = decode(enc);

        XTESTS_TEST_FAIL("should not reach here");
    }
    catch(b64::coding_exception& x)
    {
        XTESTS_TEST_ENUM_EQUAL(b64::B64_RC_DATA_ERROR, x.get_rc());
        XTESTS_TEST_POINTER_EQUAL(NULL, x.get_badChar());
    }
}

static void test_3_04(void)
{
    using namespace b64;

    stlsoft::simple_string const enc = "abcdefg%";

    try
    {
        blob_t r = decode(B64_F_STOP_ON_NOTHING, enc);

        XTESTS_TEST_PASSED();
    }
    catch(b64::coding_exception& /* x */)
    {
        XTESTS_TEST_FAIL("should not reach here");
    }
}

static void test_3_05(void)
{
    using namespace b64;

    string_t const enc = "abcdefg%";

    try
    {
        blob_t r = decode(enc);

        XTESTS_TEST_FAIL("should not reach here");
    }
    catch(b64::coding_exception& x)
    {
        XTESTS_TEST_ENUM_EQUAL(b64::B64_RC_DATA_ERROR, x.get_rc());
        XTESTS_TEST_POINTER_NOT_EQUAL(NULL, x.get_badChar());
        XTESTS_TEST_CHARACTER_EQUAL('%', *x.get_badChar());
    }
}

static void test_3_06(void)
{
    using namespace b64;

    string_t const enc = "abcdefg%";

    try
    {
        blob_t r = decode(enc, B64_F_STOP_ON_NOTHING);

        XTESTS_TEST_PASSED();
    }
    catch(b64::coding_exception& /* x */)
    {
        XTESTS_TEST_FAIL("should not reach here");
    }
}

static void test_3_07(void)
{
    using namespace b64;

    string_t const enc = "abcdefg%";

    try
    {
        blob_t r = decode(B64_F_STOP_ON_NOTHING, enc);

        XTESTS_TEST_PASSED();
    }
    catch(b64::coding_exception& /* x */)
    {
        XTESTS_TEST_FAIL("should not reach here");
    }
}

static void test_3_08(void)
{
    using namespace b64;

    string_t const enc = "abcdefg%";

    try
    {
        b64_char_t const*   badChar =   NULL;
        B64_RC              rc      =   B64_RC_OK;

        blob_t r = decode(enc, B64_F_STOP_ON_BAD_CHAR, &badChar, &rc);

        XTESTS_TEST_INTEGER_EQUAL(0u, r.size());
        XTESTS_TEST_POINTER_NOT_EQUAL(NULL, badChar);
        XTESTS_TEST_CHARACTER_EQUAL('%', *badChar);
        XTESTS_TEST_ENUM_EQUAL(B64_RC_DATA_ERROR, rc);
    }
    catch(b64::coding_exception& /* x */)
    {
        XTESTS_TEST_FAIL("should not reach here");
    }
}

static void test_3_09(void)
{
    using namespace b64;

    char const enc[] = "abcdefg%";

    try
    {
        b64_char_t const*   badChar =   NULL;
        B64_RC              rc      =   B64_RC_OK;

        blob_t r = decode(enc, STLSOFT_NUM_ELEMENTS(enc) - 1u, B64_F_STOP_ON_BAD_CHAR, &badChar, &rc);

        XTESTS_TEST_INTEGER_EQUAL(0u, r.size());
        XTESTS_TEST_POINTER_NOT_EQUAL(NULL, badChar);
        XTESTS_TEST_CHARACTER_EQUAL('%', *badChar);
        XTESTS_TEST_ENUM_EQUAL(B64_RC_DATA_ERROR, rc);
    }
    catch(b64::coding_exception& /* x */)
    {
        XTESTS_TEST_FAIL("should not reach here");
    }
}

static void test_3_10(void)
{
    using namespace b64;

    char const enc[] = "abcdefg%";

    b64_char_t const*   badChar =   NULL;

    try
    {
        blob_t r = decode(enc, STLSOFT_NUM_ELEMENTS(enc) - 1u, B64_F_STOP_ON_BAD_CHAR, &badChar);

        XTESTS_TEST_FAIL("should not reach here");
    }
    catch(b64::coding_exception& x)
    {
        XTESTS_TEST_ENUM_EQUAL(b64::B64_RC_DATA_ERROR, x.get_rc());
        XTESTS_TEST_POINTER_EQUAL(NULL, x.get_badChar());
        XTESTS_TEST_CHARACTER_EQUAL('%', *badChar);
    }
}



} /* anonymous namespace */

/* /////////////////////////////////////////////////////////////////////////
 * Test component implementation file include(s)
 */

#include <b64/implicit_link.h>

/* ///////////////////////////// end of file //////////////////////////// */
