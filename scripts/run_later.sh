echo "[$(date +'%F %T')] Sleeping for 3 hours..."
sleep 3h

echo "[$(date +'%F %T')] Starting experiment..."
python3 main.py --exp 2 quick_align
python3 main.py --exp 5 quick_align
python3 main.py --exp 6 quick_align

echo "[$(date +'%F %T')] Done."
