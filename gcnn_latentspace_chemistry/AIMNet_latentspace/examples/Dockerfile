# Use an image with Python and necessary dependencies
FROM python:3.7

# Install dependencies
RUN pip install torch==1.1.0 torchvision==0.3.0

# Set working directory
WORKDIR /app

# Copy your code into the container
COPY . /app

# Command to run your code (replace with your command)
CMD ["python", "geometry_optimization.py"]