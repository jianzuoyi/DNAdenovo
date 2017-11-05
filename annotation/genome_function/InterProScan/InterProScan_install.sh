# 1. Install InterProScan software core

wget ftp://ftp.ebi.ac.uk/pub/software/unix/iprscan/5/5.20-59.0/interproscan-5.20-59.0-64-bit.tar.gz
wget ftp://ftp.ebi.ac.uk/pub/software/unix/iprscan/5/5.20-59.0/interproscan-5.20-59.0-64-bit.tar.gz.md5

md5sum -c interproscan-5.20-59.0-64-bit.tar.gz.md5

tar -pxvzf interproscan-5.20-59.0-64-bit.tar.gz
# where:
#     p = preserve the file permissions
#     x = extract files from an archive
#     v = verbosely list the files processed
#     z = filter the archive through gzip
#     f = use archive file

# 2. Install Panther Modles

cd [InterProScan5 home]/data/

wget ftp://ftp.ebi.ac.uk/pub/software/unix/iprscan/5/data/panther-data-10.0.tar.gz
wget ftp://ftp.ebi.ac.uk/pub/software/unix/iprscan/5/data/panther-data-10.0.tar.gz.md5

md5sum -c panther-data-10.0.tar.gz.md5

tar -pxvzf panther-data-10.0.tar.gz

3. Use Pre-calculated Match Looup Server


# InterPro Homepage
http://www.ebi.ac.uk/interpro/

InterPro, integrated resource fro protein families, domains and motifs. It provides a single, consisten interface of protein signatures
contributed by ten different databases, each of wich uses a slightly different method for deriving protein signatures.






