FROM ghcr.io/astral-sh/uv:0.5-python3.10-bookworm AS build-stage

# Set environment variables
ENV UV_COMPILE_BYTECODE=1 \
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
    && uv sync --frozen --no-dev \
    && uv build

# Second stage: Runtime image
FROM python:3.10.18-bookworm AS runtime-stage

# Set environment variables
ENV BGDATA_LOCAL="/bgdatacache" \
    BBGLAB_HOME="/home/user/.config/bbglab/"

# Copy oncodrive3d from the build stage
COPY --from=build-stage /oncodrive3d/dist /oncodrive3d/dist
WORKDIR /oncodrive3d

# Install oncodrive3d
RUN pip install dist/*.tar.gz

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
RUN apt-get update && apt-get install -y curl \
    && pip install bgdata bgreference \
    && bgdata get datasets/genomereference/hg38 \
    && bgdata get datasets/genomereference/mm39 \
    && python3 -c "from bgreference import hg38; hg38(1, 1300000, 3000)" \
    && python3 -c "from bgreference import mm39; mm39(1, 1300000, 3000)" \
    && rm -rf /var/lib/apt/lists/*

# Set permissions for cache directory
RUN chmod -R 0755 "$BGDATA_LOCAL"

# Set entrypoint (optional)
CMD ["python3"]
