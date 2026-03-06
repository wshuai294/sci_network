# Deploying Science Networks

This project is a Flask app and uses port 5050 by default. You can deploy it in any of the following ways.

---

## Option 1: Run directly (quickest)

From the project root `sci_net/`:

```bash
cd /path/to/sci_net
pip install -r requirements.txt
python web/app.py
```

Open in a browser: `http://YOUR_IP:5050`. For access from other machines, keep `host="0.0.0.0"` in `app.run()` (already set) and allow port 5050 in the server firewall.

**Note:** This runs in development mode (debug=True). For production, use one of the options below.

---

## Option 2: Gunicorn + Nginx (recommended on Linux)

### 1. Install dependencies

```bash
cd /path/to/sci_net
pip install -r requirements.txt
```

### 2. Run the app with Gunicorn

Run from the **project root** (do not run from inside `web/`):

```bash
cd /path/to/sci_net
gunicorn -w 4 -b 127.0.0.1:5050 "web.wsgi:app"
```

- `-w 4`: 4 workers (adjust for your CPU)
- `-b 127.0.0.1:5050`: listen on localhost only; Nginx will reverse-proxy to it

### 3. Keep it running with systemd (optional)

Create `/etc/systemd/system/sci_net.service`:

```ini
[Unit]
Description=Science Networks Flask App
After=network.target

[Service]
User=www-data
WorkingDirectory=/path/to/sci_net
Environment="PATH=/path/to/venv/bin"
ExecStart=/path/to/venv/bin/gunicorn -w 4 -b 127.0.0.1:5050 "web.wsgi:app"
Restart=always

[Install]
WantedBy=multi-user.target
```

Then:

```bash
sudo systemctl daemon-reload
sudo systemctl enable sci_net
sudo systemctl start sci_net
```

### 4. Nginx reverse proxy

Add a server block (e.g. `/etc/nginx/sites-available/sci_net`):

```nginx
server {
    listen 80;
    server_name your-domain.com;   # change to your domain or IP

    location / {
        proxy_pass http://127.0.0.1:5050;
        proxy_set_header Host $host;
        proxy_set_header X-Real-IP $remote_addr;
        proxy_set_header X-Forwarded-For $proxy_add_x_forwarded_for;
        proxy_set_header X-Forwarded-Proto $scheme;
        proxy_read_timeout 300;
        proxy_connect_timeout 300;
        proxy_send_timeout 300;
    }
}
```

Enable and reload Nginx:

```bash
sudo ln -s /etc/nginx/sites-available/sci_net /etc/nginx/sites-enabled/
sudo nginx -t && sudo systemctl reload nginx
```

For HTTPS, use `certbot` to obtain a certificate for `your-domain.com` and configure Nginx to listen on 443.

---

## Option 3: Docker

Create a `Dockerfile` in the project root:

```dockerfile
FROM python:3.11-slim
WORKDIR /app
COPY requirements.txt .
RUN pip install --no-cache-dir -r requirements.txt
COPY . .
EXPOSE 5050
CMD ["gunicorn", "-w", "4", "-b", "0.0.0.0:5050", "web.wsgi:app"]
```

Build and run:

```bash
cd /path/to/sci_net
docker build -t sci_net .
docker run -d -p 5050:5050 --name sci_net sci_net
```

Then open `http://YOUR_SERVER_IP:5050`.

---

## Option 4: Cloud PaaS

- **Railway / Render / Fly.io**: Connect your GitHub repo, choose Python, and set the start command to:
  ```bash
  pip install -r requirements.txt && gunicorn -w 4 -b 0.0.0.0:$PORT web.wsgi:app
  ```
  Use `$PORT` if the platform provides it; otherwise use `5050`.
- **Heroku**: Create a `Procfile` in the project root:
  ```
  web: gunicorn -w 4 -b 0.0.0.0:$PORT web.wsgi:app
  ```

---

## Summary

| Option            | Use case                                      |
|-------------------|-----------------------------------------------|
| Direct run         | Local or quick demo                           |
| Gunicorn + Nginx  | Your own Linux server, domain, HTTPS          |
| Docker            | Any environment with Docker                   |
| Cloud PaaS        | No server management, use a hosted platform   |

**After deployment:** The app calls external APIs (PubMed, Google Scholar, etc.), so the server must have internet access. Large requests (e.g. 50,000 papers) can take a long time; set Nginx or platform timeouts high (e.g. 300s).
