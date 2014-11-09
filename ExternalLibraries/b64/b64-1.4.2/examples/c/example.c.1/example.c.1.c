/* /////////////////////////////////////////////////////////////////////////////
 * File:        example.c.1.c
 *
 * Purpose:     Implementation file for the example.c.1 project.
 *
 * Created:     27th April 2006
 * Updated:     21st April 2009
 *
 * Status:      Wizard-generated
 *
 * License:     (Licensed under the Synesis Software Open License)
 *
 *              Copyright (c) 2006-2009, Synesis Software Pty Ltd.
 *              All rights reserved.
 *
 *              www:        http://www.synesis.com.au/software
 *
 *              This source code is placed into the public domain 2006
 *              by Synesis Software Pty Ltd. There are no restrictions
 *              whatsoever to your use of the software. 
 *
 *              This source code is provided by Synesis Software Pty Ltd "as is"
 *              and any warranties, whether expressed or implied, including, but
 *              not limited to, the implied warranties of merchantability and
 *              fitness for a particular purpose are disclaimed. In no event
 *              shall the Synesis Software Pty Ltd be liable for any direct,
 *              indirect, incidental, special, exemplary, or consequential
 *              damages (including, but not limited to, procurement of
 *              substitute goods or services; loss of use, data, or profits; or
 *              business interruption) however caused and on any theory of
 *              liability, whether in contract, strict liability, or tort
 *              (including negligence or otherwise) arising in any way out of
 *              the use of this software, even if advised of the possibility of
 *              such damage. 
 *
 *              Neither the name of Synesis Software Pty Ltd nor the names of
 *              any subdivisions, employees or agents of Synesis Software Pty
 *              Ltd, nor the names of any other contributors to this software
 *              may be used to endorse or promote products derived from this
 *              software without specific prior written permission. 
 *
 * ////////////////////////////////////////////////////////////////////////// */


/* b64 Header Files */
#include <b64/b64.h>
#include <b64/implicit_link.h>  /* Don't include this if you want to use explicit linking */

/* Standard C Header Files */
#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

/* /////////////////////////////////////////////////////////////////////////
 * Macros
 */

#ifndef NUM_ELEMENTS
# define NUM_ELEMENTS(x)        (sizeof(x) / sizeof(0[x]))
#endif /* !NUM_ELEMENTS */

/* ////////////////////////////////////////////////////////////////////////// */

int main(void)
{
    /* Simple conversion using b64_encode() and b64_decode(). */

    /* Declare an array of bytes to use as the 'binary' blob to encode. */
    unsigned char   bytes[]  =   { 1, 2, 3, 4, 5, 6, 7, 8, 9, 10 };

    /* Invoke b64_encode() to determine the maximum size of the encoded 
     * string, specifying the argument 0 for the destSize parameter.
     * 
     * Note: The value returned may be larger than the actual size 
     * required, but will never be smaller.
     */
    size_t  cch     =   b64_encode(&bytes[0], sizeof(bytes), NULL, 0);

    /* Using the length determined by the call to b64_encode(), create 
     * a buffer of sufficient size.
     */
    char    *enc    =   (char*)malloc(cch);
    size_t  cb;
    void    *dec;
    size_t  i;

    printf("Converting %u bytes:\n", (unsigned)NUM_ELEMENTS(bytes));
    for(i = 0; i != NUM_ELEMENTS(bytes); ++i)
    {
        printf(" %d", bytes[i]);
    }
    puts("");

    /* Perform the encoding. */
    cch = b64_encode(&bytes[0], sizeof(bytes), enc, cch);

    printf("Encoded form: [%.*s]\n", (int)cch, enc);

    /* Invoke b64_decode() to determine the maximum size of the decoded
     * blob, specifying the argument 0 for the destLen parameter.
     *
     * Note: The value returned may be larger than the actual size
     * required, but will never be smaller.
     */
    cb = b64_decode(enc, cch, NULL, 0);

    /* Using the length determined by the call to b64_decode(), create 
     * a buffer of sufficient size.
     */
    dec = malloc(cb);

    /* Perform the decoding. */
    cb = b64_decode(enc, cch, dec, cb);

    /* Verify that the decoding is exactly the same size and contents as 
     * the encoding.
     */
    assert(cb == sizeof(bytes));
    assert(0 == memcmp(&bytes[0], dec, sizeof(bytes)));

    /* Clean up */
    free(dec);
    free(enc);

    /* Placate the compiler's warnings. */
    ((void)cb);

    return EXIT_SUCCESS;
}


/* ////////////////////////////////////////////////////////////////////////// */
