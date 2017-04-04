# <p >GSE 2.x or IMS 1.0 INT and CM6 waveform data to miniSEED converter</p>

1. [Name](#)
1. [Synopsis](#synopsis)
1. [Description](#description)
1. [Options](#options)
1. [List Files](#list-files)
1. [About Gse & Ims Formats](#about-gse-&-ims-formats)
1. [Author](#author)

## <a id='synopsis'>Synopsis</a>

<pre >
gse2mseed [options] file1 [file2 file3 ...]
</pre>

## <a id='description'>Description</a>

<p ><b>gse2mseed</b> converts GSE 2.x or IMS 1.0 INT or CM6 compressed waveform data files to miniSEED.  One or more input files may be specified on the command line.  If an input file name is prefixed with an '@' character the file is assumed to contain a list of input data files, see <i>LIST FILES</i> below.</p>

<p >By default all data from a given input file is written to a file of the same name with a ".mseed" suffix.  If the input file name includes a ".gse" suffix/extension it will be removed.  The output data may be re-directed to a single file or stdout using the -o option.</p>

## <a id='options'>Options</a>

<b>-V</b>

<p style="padding-left: 30px;">Print program version and exit.</p>

<b>-h</b>

<p style="padding-left: 30px;">Print program usage and exit.</p>

<b>-v</b>

<p style="padding-left: 30px;">Be more verbose.  This flag can be used multiple times ("-v -v" or "-vv") for more verbosity.</p>

<b>-i</b>

<p style="padding-left: 30px;">Ignore GSE checksum mismatch errors, a warning message is printed but processing continues.</p>

<b>-B</b>

<p style="padding-left: 30px;">Buffer all input data into memory before packing it into miniSEED records.  The host computer must have enough memory to store all of the data.  By default the program will flush it's data buffers after each input block is read.  An output file must be specified with the http://www.seismo.ethz.ch/autodrm/ -o option when using this option.</p>

<b>-n </b><i>netcode</i>

<p style="padding-left: 30px;">Specify the SEED network code to use, if not specified the network code will be blank.  It is highly recommended to specify a network code.</p>

<b>-l </b><i>loccode</i>

<p style="padding-left: 30px;">Specify the SEED location ID to use, if not specified the GSE auxiluarly idenficiation code will be used.</p>

<b>-r </b><i>bytes</i>

<p style="padding-left: 30px;">Specify the miniSEED record length in <i>bytes</i>, default is 4096.</p>

<b>-e </b><i>encoding</i>

<p style="padding-left: 30px;">Specify the miniSEED data encoding format, default is 11 (Steim-2 compression).  Other supported encoding formats include 10 (Steim-1 compression), 1 (16-bit integers) and 3 (32-bit integers).  The 16-bit integers encoding should only be used if all data samples can be represented in 16 bits.</p>

<b>-b </b><i>byteorder</i>

<p style="padding-left: 30px;">Specify the miniSEED byte order, default is 1 (big-endian or most significant byte first).  The other option is 0 (little-endian or least significant byte first).  It is highly recommended to always create big-endian SEED.</p>

<b>-o </b><i>outfile</i>

<p style="padding-left: 30px;">Write all miniSEED records to <i>outfile</i>, if <i>outfile</i> is a single dash (-) then all miniSEED output will go to stdout.  All diagnostic output from the program is written to stderr and should never get mixed with data going to stdout.</p>

## <a id='list-files'>List Files</a>

<p >If an input file is prefixed with an '@' character the file is assumed to contain a list of file for input.  The list should be a simple text file with one input file name per line.</p>

<p >Multiple list files can be combined with multiple input files on the command line.</p>

<p >An example of a simple test list:</p>

<pre >
SENIN.CH.gse
ANMO.IU.gse
</pre>

## <a id='about-gse-&-ims-formats'>About Gse & Ims Formats</a>

<p >The GSE format was created by the Group of Scientific Experts during a series of Technical Tests.  With minor modification the GSE 2.x format was adopted by the Comprehensive Nuclear-Test-Ban Treaty Organization's (CTBTO) International Monitoring System and designated IMS 1.0.</p>

<p >The GSE format is commonly used in the Swiss AutoDRM: http://www.seismo.ethz.ch/prod/autodrm/index_EN</p>

## <a id='author'>Author</a>

<pre >
Chad Trabant
IRIS Data Management Center
</pre>


(man page 2017/04/03)
