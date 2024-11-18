FROM python:3.8

# Set environment variables
ENV LC_ALL=C \
    BGDATA_LOCAL="/bgdatacache" \
    BBGLAB_HOME="/home/user/.config/bbglab/"

# Update system and install necessary tools
RUN apt-get update -y \
    && apt-get upgrade -y \
    && apt-get install -y --no-install-recommends bash \
    && pip install --upgrade pip \
    && apt-get clean \
    && rm -rf /var/lib/apt/lists/*

# Create necessary directories
RUN mkdir -p "$BGDATA_LOCAL" "$BBGLAB_HOME"

# Install bgdata and bgreference packages
RUN pip install bgdata bgreference 

# Write bgdata configuration
RUN echo "# Version of the bgdata config file\n\
version = 2\n\
\n\
# The default local folder to store the data packages\n\
local_repository = \"$BGDATA_LOCAL\"\n\
\n\
# The remote URL to download the data packages\n\
remote_repository = \"https://bbglab.irbbarcelona.org/bgdata\"\n\
\n\
# If you want to force bgdata to work only locally\n\
# offline = True\n\
\n\
# Cache repositories\n\
[cache_repositories]" > "$BBGLAB_HOME/bgdatav2.conf"

# Pre-fetch and prepare genome data
RUN bgdata get datasets/genomereference/hg38 && bgdata get datasets/genomereference/hg19 \
    && bgdata get datasets/genomereference/mm39 && bgdata get datasets/genomereference/mm10 \
    && python3 -c "from bgreference import hg38, hg19; hg38(1, 1300000, 3000); hg19(1, 1300000, 3000)" \
    && python3 -c "from bgreference import mm39, mm10; mm39(1, 1300000, 3000); mm10(1, 1300000, 3000)"

# Adjust permissions for bgdatacache
RUN chmod -R 777 "$BGDATA_LOCAL"

# Copy and install Oncodrive3D
COPY oncodrive3d /oncodrive3d
RUN pip install /oncodrive3d

# Metadata
LABEL Author="stefano.pellegrini@irbbarcelona.org"

# Test the installation
CMD ["oncodrive3D", "run", "--help"]