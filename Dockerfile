# Use the uv image with Python 3.10 based on Bookworm slim
FROM ghcr.io/astral-sh/uv:0.5-python3.10-bookworm

# Set environment variables
ENV UV_COMPILE_BYTECODE=1 \
    BGDATA_LOCAL="/bgdatacache" \
    BBGLAB_HOME="/home/user/.config/bbglab/"
 
# Set the working directory to /oncodrive3d
WORKDIR /oncodrive3d

# Stage necessary files into the container
COPY uv.lock uv.lock
COPY pyproject.toml pyproject.toml
COPY scripts scripts
COPY README.md README.md

# Install dependencies and build the project
RUN mkdir -p /root/.cache/uv \
    && cd /oncodrive3d \
    && uv sync --frozen --no-dev \
    && uv build \
    && pip install dist/*.tar.gz \
    # Clean up unnecessary files
    && rm -rf /root/.cache/uv /oncodrive3d/build /oncodrive3d/.venv

# Create required directories
RUN install -d -m 0755 "$BGDATA_LOCAL" "$BBGLAB_HOME"
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

# Set permissions for cache directory
RUN chmod -R 0755 "$BGDATA_LOCAL"