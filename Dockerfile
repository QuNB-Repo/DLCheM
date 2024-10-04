# Use a lightweight Python image as the base
FROM python:3.7-slim

# Set the working directory in the container
WORKDIR /app

# Copy the requirements file to install any dependencies
COPY requirements.txt /app/

# Install dependencies
RUN pip install --no-cache-dir -r requirements.txt

# Copy the rest of your application's code
COPY . /app

# Expose the port that the application runs on (adjust this if needed)
EXPOSE 5000

# Command to run your application (adjust to your entrypoint script)
CMD ["python", "app.py"]

