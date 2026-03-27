#!/usr/bin/env bash
# Deploy text2gene2 to a fresh Debian/Ubuntu VPS.
# Run as root. Tested on Ubuntu 24.04 LTS.
set -euo pipefail

APP_DIR=/opt/text2gene2
REPO_URL=https://github.com/nthmost/text2gene2.git   # update if needed

echo "==> Installing system packages..."
apt-get update -qq
apt-get install -y python3.12 python3.12-venv python3-pip redis-server nginx certbot python3-certbot-nginx git

echo "==> Enabling and starting Redis..."
systemctl enable redis-server
systemctl start redis-server

echo "==> Cloning / updating repo..."
if [ -d "$APP_DIR/.git" ]; then
    git -C "$APP_DIR" pull
else
    git clone "$REPO_URL" "$APP_DIR"
fi

echo "==> Setting up Python venv..."
python3.12 -m venv "$APP_DIR/venv"
"$APP_DIR/venv/bin/pip" install --quiet --upgrade pip
"$APP_DIR/venv/bin/pip" install --quiet "$APP_DIR"

echo "==> Installing systemd service..."
cp "$APP_DIR/deploy/text2gene2.service" /etc/systemd/system/
systemctl daemon-reload
systemctl enable text2gene2

echo "==> Installing nginx config..."
cp "$APP_DIR/deploy/nginx.conf" /etc/nginx/sites-available/text2gene2
ln -sf /etc/nginx/sites-available/text2gene2 /etc/nginx/sites-enabled/text2gene2
rm -f /etc/nginx/sites-enabled/default
nginx -t && systemctl reload nginx

echo ""
echo "==> DONE. Next steps:"
echo "    1. Copy .env.example to $APP_DIR/.env and fill in API keys"
echo "    2. chown -R www-data:www-data $APP_DIR"
echo "    3. systemctl start text2gene2"
echo "    4. certbot --nginx -d text2gene.org -d www.text2gene.org"
