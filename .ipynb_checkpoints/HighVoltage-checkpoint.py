{
 "cells": [
  {
   "cell_type": "code",
   "execution_count": 2,
   "metadata": {},
   "outputs": [],
   "source": [
    "import overview\n",
    "import importlib\n",
    "import rawFileReader\n",
    "import pickle\n",
    "from matplotlib.backends.backend_pdf import PdfPages\n",
    "import pandas as pd\n",
    "\n",
    "dir_path = \"/usera/jd2052/Documents/anubis/\" # insert your directory path\n"
   ]
  },
  {
   "cell_type": "code",
   "execution_count": 7,
   "metadata": {},
   "outputs": [
    {
     "name": "stdout",
     "output_type": "stream",
     "text": [
      "All good\n"
     ]
    }
   ],
   "source": [
    "#Save chunks\n",
    "files = [\"proAnubis_240815_1759.raw\", \"proAnubis_240818_1125.raw\",\"proAnubis_240818_1125.raw\", \"proAnubis_240818_1325.raw\", \"proAnubis_240829_2026.raw\"]\n",
    "rpc = 4\n",
    "\n",
    "storage_name = f\"chunks_hv{rpc}.pkl\"\n",
    "file_name = files[rpc]\n",
    "\n",
    "hv_file = 'data/hvScan.csv'  # Replace with your file path\n",
    "hv_data = pd.read_csv(hv_file, usecols=[\"start_\"+str(rpc),\"end_\"+str(rpc),\"voltage_\"+str(rpc)])\n",
    "#rpc 0: 2024-08-15 17:24:22\t\n",
    "#rpc 1\t2024-08-18 11:39:42\n",
    "#rpc 2\t2024-08-18 12:22:21 \n"
   ]
  },
  {
   "cell_type": "code",
   "execution_count": 9,
   "metadata": {},
   "outputs": [
    {
     "name": "stdout",
     "output_type": "stream",
     "text": [
      "Voltage: 4000\n",
      "Initial time 2024-08-29 19:26:05.663040 1724955965.66304\n",
      "Start: 1724957390.0\n",
      "Bigger than start: False\n"
     ]
    },
    {
     "name": "stderr",
     "output_type": "stream",
     "text": [
      "Skipping Events proAnubis_240829_2026.raw: 1610Events [20:36,  1.30Events/s]    \n"
     ]
    },
    {
     "name": "stdout",
     "output_type": "stream",
     "text": [
      "Skip time 2024-08-29 19:53:02.892229\n"
     ]
    },
    {
     "name": "stderr",
     "output_type": "stream",
     "text": [
      "Processing Chunks proAnubis_240829_2026.raw: 197.21seconds [00:00, 2987.83seconds/s]\n"
     ]
    },
    {
     "name": "stdout",
     "output_type": "stream",
     "text": [
      "Ending time: 2024-08-29 19:53:07.213203\n",
      "Voltage: 4200\n",
      "Initial time 2024-08-29 19:53:13.705631 1724957593.705631\n",
      "Start: 1724957437.0\n",
      "Bigger than start: False\n"
     ]
    },
    {
     "name": "stderr",
     "output_type": "stream",
     "text": [
      "Skipping Events proAnubis_240829_2026.raw: 0Events [00:00, ?Events/s]\n"
     ]
    },
    {
     "name": "stdout",
     "output_type": "stream",
     "text": [
      "Skip time 2024-08-29 19:53:25.029169\n"
     ]
    },
    {
     "name": "stderr",
     "output_type": "stream",
     "text": [
      "Processing Chunks proAnubis_240829_2026.raw: 173.35seconds [00:00, 2766.32seconds/s]\n"
     ]
    },
    {
     "name": "stdout",
     "output_type": "stream",
     "text": [
      "Ending time: 2024-08-29 19:53:30.353207\n",
      "Voltage: 4300\n",
      "Initial time 2024-08-29 19:53:36.494910 1724957616.49491\n",
      "Start: 1724957488.0\n",
      "Bigger than start: False\n"
     ]
    },
    {
     "name": "stderr",
     "output_type": "stream",
     "text": [
      "Skipping Events proAnubis_240829_2026.raw: 0Events [00:00, ?Events/s]\n"
     ]
    },
    {
     "name": "stdout",
     "output_type": "stream",
     "text": [
      "Skip time 2024-08-29 19:53:41.576758\n"
     ]
    },
    {
     "name": "stderr",
     "output_type": "stream",
     "text": [
      "Processing Chunks proAnubis_240829_2026.raw: 137.36seconds [00:00, 2910.42seconds/s]\n"
     ]
    },
    {
     "name": "stdout",
     "output_type": "stream",
     "text": [
      "Ending time: 2024-08-29 19:53:45.363551\n",
      "Voltage: 4400\n",
      "Initial time 2024-08-29 19:53:49.942240 1724957629.94224\n",
      "Start: 1724957535.0\n",
      "Bigger than start: False\n"
     ]
    },
    {
     "name": "stderr",
     "output_type": "stream",
     "text": [
      "Skipping Events proAnubis_240829_2026.raw: 0Events [00:00, ?Events/s]\n"
     ]
    },
    {
     "name": "stdout",
     "output_type": "stream",
     "text": [
      "Skip time 2024-08-29 19:53:53.140415\n"
     ]
    },
    {
     "name": "stderr",
     "output_type": "stream",
     "text": [
      "Processing Chunks proAnubis_240829_2026.raw: 100.46seconds [00:00, 1694.07seconds/s]"
     ]
    },
    {
     "name": "stdout",
     "output_type": "stream",
     "text": [
      "Ending time: 2024-08-29 19:53:55.459334\n",
      "Voltage: 4500\n"
     ]
    },
    {
     "name": "stderr",
     "output_type": "stream",
     "text": [
      "\n"
     ]
    },
    {
     "name": "stdout",
     "output_type": "stream",
     "text": [
      "Initial time 2024-08-29 19:53:59.803030 1724957639.80303\n",
      "Start: 1724957573.0\n",
      "Bigger than start: False\n"
     ]
    },
    {
     "name": "stderr",
     "output_type": "stream",
     "text": [
      "Skipping Events proAnubis_240829_2026.raw: 0Events [00:00, ?Events/s]\n"
     ]
    },
    {
     "name": "stdout",
     "output_type": "stream",
     "text": [
      "Skip time 2024-08-29 19:54:04.869307\n"
     ]
    },
    {
     "name": "stderr",
     "output_type": "stream",
     "text": [
      "Processing Chunks proAnubis_240829_2026.raw: 52.58seconds [00:00, 1236.39seconds/s]\n"
     ]
    },
    {
     "name": "stdout",
     "output_type": "stream",
     "text": [
      "Ending time: 2024-08-29 19:53:45.575269\n",
      "Voltage: 4600\n",
      "Initial time 2024-08-29 19:54:12.567272 1724957652.567272\n",
      "Start: 1724957612.0\n",
      "Bigger than start: False\n"
     ]
    },
    {
     "name": "stderr",
     "output_type": "stream",
     "text": [
      "Skipping Events proAnubis_240829_2026.raw: 0Events [00:00, ?Events/s]\n"
     ]
    },
    {
     "name": "stdout",
     "output_type": "stream",
     "text": [
      "Skip time 2024-08-29 19:54:15.146708\n"
     ]
    },
    {
     "name": "stderr",
     "output_type": "stream",
     "text": [
      "Processing Chunks proAnubis_240829_2026.raw: 46.1seconds [00:00, 1284.56seconds/s]\n"
     ]
    },
    {
     "name": "stdout",
     "output_type": "stream",
     "text": [
      "Ending time: 2024-08-29 19:54:18.103638\n",
      "Voltage: 4700\n",
      "Initial time 2024-08-29 19:54:21.008595 1724957661.008595\n",
      "Start: 1724957653.0\n",
      "Bigger than start: False\n"
     ]
    },
    {
     "name": "stderr",
     "output_type": "stream",
     "text": [
      "Skipping Events proAnubis_240829_2026.raw: 0Events [00:00, ?Events/s]\n"
     ]
    },
    {
     "name": "stdout",
     "output_type": "stream",
     "text": [
      "Skip time 2024-08-29 19:54:22.337390\n"
     ]
    },
    {
     "name": "stderr",
     "output_type": "stream",
     "text": [
      "Processing Chunks proAnubis_240829_2026.raw: 33.50999999999999seconds [00:00, 77.54seconds/s]\n"
     ]
    },
    {
     "name": "stdout",
     "output_type": "stream",
     "text": [
      "Ending time: 2024-08-29 19:54:46.511968\n",
      "Voltage: 4800\n",
      "Initial time 2024-08-29 19:54:48.918308 1724957688.918308\n",
      "Start: 1724957696.0\n",
      "Bigger than start: False\n"
     ]
    },
    {
     "name": "stderr",
     "output_type": "stream",
     "text": [
      "Skipping Events proAnubis_240829_2026.raw: 26Events [00:00, 79.90Events/s]      \n"
     ]
    },
    {
     "name": "stdout",
     "output_type": "stream",
     "text": [
      "Skip time 2024-08-29 19:55:15.359006\n"
     ]
    },
    {
     "name": "stderr",
     "output_type": "stream",
     "text": [
      "Processing Chunks proAnubis_240829_2026.raw: 72.64999999999999seconds [00:01, 65.95seconds/s]\n"
     ]
    },
    {
     "name": "stdout",
     "output_type": "stream",
     "text": [
      "Ending time: 2024-08-29 19:56:08.652910\n",
      "Voltage: 4900\n",
      "Initial time 2024-08-29 19:56:09.474734 1724957769.474734\n",
      "Start: 1724957779.0\n",
      "Bigger than start: False\n"
     ]
    },
    {
     "name": "stderr",
     "output_type": "stream",
     "text": [
      "Skipping Events proAnubis_240829_2026.raw: 11Events [00:00, 29.06Events/s]      \n"
     ]
    },
    {
     "name": "stdout",
     "output_type": "stream",
     "text": [
      "Skip time 2024-08-29 19:56:20.989487\n"
     ]
    },
    {
     "name": "stderr",
     "output_type": "stream",
     "text": [
      "Processing Chunks proAnubis_240829_2026.raw:  96%|▉| 57.739999999999995/60.0 [00/usera/jd2052/.local/lib/python3.9/site-packages/tqdm/std.py:636: TqdmWarning: clamping frac to range [0, 1]\n",
      "  full_bar = Bar(frac,\n",
      "Processing Chunks proAnubis_240829_2026.raw: 100%|█| 60.13999999999999/60.0 [00:\n"
     ]
    },
    {
     "name": "stdout",
     "output_type": "stream",
     "text": [
      "Ending time: 2024-08-29 19:57:19.142044\n",
      "Voltage: 5000\n",
      "Initial time 2024-08-29 19:57:19.723654 1724957839.723654\n",
      "Start: 1724957850.0\n",
      "Bigger than start: False\n"
     ]
    },
    {
     "name": "stderr",
     "output_type": "stream",
     "text": [
      "Skipping Events proAnubis_240829_2026.raw: 11Events [00:00, 15.00Events/s]      \n"
     ]
    },
    {
     "name": "stdout",
     "output_type": "stream",
     "text": [
      "Skip time 2024-08-29 19:57:31.272498\n"
     ]
    },
    {
     "name": "stderr",
     "output_type": "stream",
     "text": [
      "Processing Chunks proAnubis_240829_2026.raw: 100%|█| 31.119999999999976/31.0 [00\n"
     ]
    },
    {
     "name": "stdout",
     "output_type": "stream",
     "text": [
      "Ending time: 2024-08-29 19:58:01.119955\n",
      "Voltage: 5100\n",
      "Initial time 2024-08-29 19:58:01.239046 1724957881.239046\n",
      "Start: 1724957891.0\n",
      "Bigger than start: False\n"
     ]
    },
    {
     "name": "stderr",
     "output_type": "stream",
     "text": [
      "Skipping Events proAnubis_240829_2026.raw: 100%|█| 10/10 [00:01<00:00,  7.92Even\n"
     ]
    },
    {
     "name": "stdout",
     "output_type": "stream",
     "text": [
      "Skip time 2024-08-29 19:58:11.299222\n"
     ]
    },
    {
     "name": "stderr",
     "output_type": "stream",
     "text": [
      "Processing Chunks proAnubis_240829_2026.raw: 100%|█| 32.009999999999984/32.0 [00\n"
     ]
    },
    {
     "name": "stdout",
     "output_type": "stream",
     "text": [
      "Ending time: 2024-08-29 19:58:43.009870\n",
      "Voltage: 5200\n",
      "Initial time 2024-08-29 19:58:43.379558 1724957923.379558\n",
      "Start: 1724957931.0\n",
      "Bigger than start: False\n"
     ]
    },
    {
     "name": "stderr",
     "output_type": "stream",
     "text": [
      "Skipping Events proAnubis_240829_2026.raw: 100%|█| 8/8 [00:01<00:00,  7.31Events\n"
     ]
    },
    {
     "name": "stdout",
     "output_type": "stream",
     "text": [
      "Skip time 2024-08-29 19:58:51.465933\n"
     ]
    },
    {
     "name": "stderr",
     "output_type": "stream",
     "text": [
      "Processing Chunks proAnubis_240829_2026.raw: 100%|█| 30.029999999999973/30.0 [00\n"
     ]
    },
    {
     "name": "stdout",
     "output_type": "stream",
     "text": [
      "Ending time: 2024-08-29 19:59:21.031913\n",
      "Voltage: 5300\n",
      "Initial time 2024-08-29 19:59:13.758136 1724957953.758136\n",
      "Start: 1724957969.0\n",
      "Bigger than start: False\n"
     ]
    },
    {
     "name": "stderr",
     "output_type": "stream",
     "text": [
      "Skipping Events proAnubis_240829_2026.raw: 100%|█| 15/15 [00:01<00:00, 10.67Even\n"
     ]
    },
    {
     "name": "stdout",
     "output_type": "stream",
     "text": [
      "Skip time 2024-08-29 19:59:29.279952\n"
     ]
    },
    {
     "name": "stderr",
     "output_type": "stream",
     "text": [
      "Processing Chunks proAnubis_240829_2026.raw: 100%|█| 31.05999999999998/31.0 [00:\n"
     ]
    },
    {
     "name": "stdout",
     "output_type": "stream",
     "text": [
      "Ending time: 2024-08-29 20:00:00.058570\n",
      "Voltage: 5400\n",
      "Initial time 2024-08-29 20:00:00.116024 1724958000.116024\n",
      "Start: 1724958009.0\n",
      "Bigger than start: False\n"
     ]
    },
    {
     "name": "stderr",
     "output_type": "stream",
     "text": [
      "Skipping Events proAnubis_240829_2026.raw: 100%|█| 9/9 [00:01<00:00,  5.03Events\n"
     ]
    },
    {
     "name": "stdout",
     "output_type": "stream",
     "text": [
      "Skip time 2024-08-29 20:00:09.658058\n"
     ]
    },
    {
     "name": "stderr",
     "output_type": "stream",
     "text": [
      "Processing Chunks proAnubis_240829_2026.raw: 100%|█| 30.109999999999996/30.0 [00\n"
     ]
    },
    {
     "name": "stdout",
     "output_type": "stream",
     "text": [
      "Ending time: 2024-08-29 20:00:39.110326\n",
      "Voltage: 5450\n",
      "Initial time 2024-08-29 20:00:39.255391 1724958039.255391\n",
      "Start: 1724958048.0\n",
      "Bigger than start: False\n"
     ]
    },
    {
     "name": "stderr",
     "output_type": "stream",
     "text": [
      "Skipping Events proAnubis_240829_2026.raw: 100%|█| 9/9 [00:04<00:00,  1.82Events\n"
     ]
    },
    {
     "name": "stdout",
     "output_type": "stream",
     "text": [
      "Skip time 2024-08-29 20:00:48.420546\n"
     ]
    },
    {
     "name": "stderr",
     "output_type": "stream",
     "text": [
      "Processing Chunks proAnubis_240829_2026.raw: 100%|█| 40.05999999999998/40.0 [00:\n"
     ]
    },
    {
     "name": "stdout",
     "output_type": "stream",
     "text": [
      "Ending time: 2024-08-29 20:01:28.056450\n",
      "Voltage: 5500\n",
      "Initial time 2024-08-29 20:01:28.116428 1724958088.116428\n",
      "Start: 1724958095.0\n",
      "Bigger than start: False\n"
     ]
    },
    {
     "name": "stderr",
     "output_type": "stream",
     "text": [
      "Skipping Events proAnubis_240829_2026.raw: 100%|█| 7/7 [00:01<00:00,  4.30Events\n"
     ]
    },
    {
     "name": "stdout",
     "output_type": "stream",
     "text": [
      "Skip time 2024-08-29 20:01:35.346569\n"
     ]
    },
    {
     "name": "stderr",
     "output_type": "stream",
     "text": [
      "Processing Chunks proAnubis_240829_2026.raw: 101%|█| 31.349999999999987/31.0 [00\n"
     ]
    },
    {
     "name": "stdout",
     "output_type": "stream",
     "text": [
      "Ending time: 2024-08-29 20:02:06.345947\n",
      "Voltage: 5550\n",
      "Initial time 2024-08-29 20:02:06.385995 1724958126.385995\n",
      "Start: 1724958135.0\n",
      "Bigger than start: False\n"
     ]
    },
    {
     "name": "stderr",
     "output_type": "stream",
     "text": [
      "Skipping Events proAnubis_240829_2026.raw: 100%|█| 9/9 [00:01<00:00,  5.02Events\n"
     ]
    },
    {
     "name": "stdout",
     "output_type": "stream",
     "text": [
      "Skip time 2024-08-29 20:02:15.233586\n"
     ]
    },
    {
     "name": "stderr",
     "output_type": "stream",
     "text": [
      "Processing Chunks proAnubis_240829_2026.raw: 100%|█| 31.030000000000005/31.0 [00\n"
     ]
    },
    {
     "name": "stdout",
     "output_type": "stream",
     "text": [
      "Ending time: 2024-08-29 20:02:46.028845\n",
      "Voltage: 5600\n",
      "Initial time 2024-08-29 20:02:46.081632 1724958166.081632\n",
      "Start: 1724958173.0\n",
      "Bigger than start: False\n"
     ]
    },
    {
     "name": "stderr",
     "output_type": "stream",
     "text": [
      "Skipping Events proAnubis_240829_2026.raw: 100%|█| 7/7 [00:01<00:00,  4.22Events\n"
     ]
    },
    {
     "name": "stdout",
     "output_type": "stream",
     "text": [
      "Skip time 2024-08-29 20:02:53.345828\n"
     ]
    },
    {
     "name": "stderr",
     "output_type": "stream",
     "text": [
      "Processing Chunks proAnubis_240829_2026.raw: 100%|█| 32.08999999999996/32.0 [00:\n"
     ]
    },
    {
     "name": "stdout",
     "output_type": "stream",
     "text": [
      "Ending time: 2024-08-29 20:03:25.085782\n",
      "Voltage: 5650\n",
      "Initial time 2024-08-29 20:03:25.108050 1724958205.10805\n",
      "Start: 1724958212.0\n",
      "Bigger than start: False\n"
     ]
    },
    {
     "name": "stderr",
     "output_type": "stream",
     "text": [
      "Skipping Events proAnubis_240829_2026.raw: 100%|█| 7/7 [00:02<00:00,  3.25Events\n"
     ]
    },
    {
     "name": "stdout",
     "output_type": "stream",
     "text": [
      "Skip time 2024-08-29 20:03:32.574447\n"
     ]
    },
    {
     "name": "stderr",
     "output_type": "stream",
     "text": [
      "Processing Chunks proAnubis_240829_2026.raw: 100%|█| 31.02999999999998/31.0 [00:\n"
     ]
    },
    {
     "name": "stdout",
     "output_type": "stream",
     "text": [
      "Ending time: 2024-08-29 20:04:03.032812\n",
      "Voltage: 5700\n",
      "Initial time 2024-08-29 20:04:03.405521 1724958243.405521\n",
      "Start: 1724958251.0\n",
      "Bigger than start: False\n"
     ]
    },
    {
     "name": "stderr",
     "output_type": "stream",
     "text": [
      "Skipping Events proAnubis_240829_2026.raw: 100%|█| 8/8 [00:02<00:00,  3.32Events\n"
     ]
    },
    {
     "name": "stdout",
     "output_type": "stream",
     "text": [
      "Skip time 2024-08-29 20:04:11.460635\n"
     ]
    },
    {
     "name": "stderr",
     "output_type": "stream",
     "text": [
      "Processing Chunks proAnubis_240829_2026.raw: 100%|█| 31.009999999999984/31.0 [00\n"
     ]
    },
    {
     "name": "stdout",
     "output_type": "stream",
     "text": [
      "Ending time: 2024-08-29 20:04:42.009502\n",
      "Voltage: 5750\n",
      "Initial time 2024-08-29 20:04:42.046315 1724958282.046315\n",
      "Start: 1724958289.0\n",
      "Bigger than start: False\n"
     ]
    },
    {
     "name": "stderr",
     "output_type": "stream",
     "text": [
      "Skipping Events proAnubis_240829_2026.raw: 8Events [00:02,  2.69Events/s]       \n"
     ]
    },
    {
     "name": "stdout",
     "output_type": "stream",
     "text": [
      "Skip time 2024-08-29 20:04:50.533769\n"
     ]
    },
    {
     "name": "stderr",
     "output_type": "stream",
     "text": [
      "Processing Chunks proAnubis_240829_2026.raw: 100%|█| 41.199999999999974/41.0 [00\n"
     ]
    },
    {
     "name": "stdout",
     "output_type": "stream",
     "text": [
      "Ending time: 2024-08-29 20:05:30.201777\n",
      "Voltage: 5800\n",
      "Initial time 2024-08-29 20:05:30.307997 1724958330.307997\n",
      "Start: 1724958341.0\n",
      "Bigger than start: False\n"
     ]
    },
    {
     "name": "stderr",
     "output_type": "stream",
     "text": [
      "Skipping Events proAnubis_240829_2026.raw: 100%|█| 11/11 [00:03<00:00,  3.44Even\n"
     ]
    },
    {
     "name": "stdout",
     "output_type": "stream",
     "text": [
      "Skip time 2024-08-29 20:05:41.512626\n"
     ]
    },
    {
     "name": "stderr",
     "output_type": "stream",
     "text": [
      "Processing Chunks proAnubis_240829_2026.raw: 101%|█| 30.25999999999999/30.0 [00:\n"
     ]
    },
    {
     "name": "stdout",
     "output_type": "stream",
     "text": [
      "Ending time: 2024-08-29 20:06:11.257813\n",
      "Voltage: 5850\n",
      "Initial time 2024-08-29 20:06:11.324510 1724958371.32451\n",
      "Start: 1724958381.0\n",
      "Bigger than start: False\n"
     ]
    },
    {
     "name": "stderr",
     "output_type": "stream",
     "text": [
      "Skipping Events proAnubis_240829_2026.raw: 100%|█| 10/10 [00:02<00:00,  3.74Even\n"
     ]
    },
    {
     "name": "stdout",
     "output_type": "stream",
     "text": [
      "Skip time 2024-08-29 20:06:21.096317\n"
     ]
    },
    {
     "name": "stderr",
     "output_type": "stream",
     "text": [
      "Processing Chunks proAnubis_240829_2026.raw: 100%|█| 31.099999999999987/31.0 [00\n"
     ]
    },
    {
     "name": "stdout",
     "output_type": "stream",
     "text": [
      "Ending time: 2024-08-29 20:06:52.104017\n",
      "Voltage: 5900\n",
      "Initial time 2024-08-29 20:06:52.294613 1724958412.294613\n",
      "Start: 1724958421.0\n",
      "Bigger than start: False\n"
     ]
    },
    {
     "name": "stderr",
     "output_type": "stream",
     "text": [
      "Skipping Events proAnubis_240829_2026.raw: 100%|█| 9/9 [00:03<00:00,  2.59Events\n"
     ]
    },
    {
     "name": "stdout",
     "output_type": "stream",
     "text": [
      "Skip time 2024-08-29 20:07:01.875485\n"
     ]
    },
    {
     "name": "stderr",
     "output_type": "stream",
     "text": [
      "Processing Chunks proAnubis_240829_2026.raw: 100%|█| 32.07/32.0 [00:25<00:00,  1\n"
     ]
    },
    {
     "name": "stdout",
     "output_type": "stream",
     "text": [
      "Ending time: 2024-08-29 20:07:33.073927\n",
      "Voltage: 5950\n",
      "Initial time 2024-08-29 20:07:33.518146 1724958453.518146\n",
      "Start: 1724958484.0\n",
      "Bigger than start: False\n"
     ]
    },
    {
     "name": "stderr",
     "output_type": "stream",
     "text": [
      "Skipping Events proAnubis_240829_2026.raw: 31Events [00:25,  1.20Events/s]      \n"
     ]
    },
    {
     "name": "stdout",
     "output_type": "stream",
     "text": [
      "Skip time 2024-08-29 20:08:04.501970\n"
     ]
    },
    {
     "name": "stderr",
     "output_type": "stream",
     "text": [
      "Processing Chunks proAnubis_240829_2026.raw: 100%|█| 32.07/32.0 [00:23<00:00,  1\n"
     ]
    },
    {
     "name": "stdout",
     "output_type": "stream",
     "text": [
      "Ending time: 2024-08-29 20:08:36.065489\n",
      "Voltage: 6000\n",
      "Initial time 2024-08-29 20:08:36.391279 1724958516.391279\n",
      "Start: 1724958537.0\n",
      "Bigger than start: False\n"
     ]
    },
    {
     "name": "stderr",
     "output_type": "stream",
     "text": [
      "Skipping Events proAnubis_240829_2026.raw: 100%|█| 21/21 [00:23<00:00,  1.10s/Ev\n"
     ]
    },
    {
     "name": "stdout",
     "output_type": "stream",
     "text": [
      "Skip time 2024-08-29 20:08:57.563492\n"
     ]
    },
    {
     "name": "stderr",
     "output_type": "stream",
     "text": [
      "Processing Chunks proAnubis_240829_2026.raw: 100%|█| 32.129999999999974/32.0 [00"
     ]
    },
    {
     "name": "stdout",
     "output_type": "stream",
     "text": [
      "Ending time: 2024-08-29 20:09:29.125644\n"
     ]
    },
    {
     "name": "stderr",
     "output_type": "stream",
     "text": [
      "\n"
     ]
    }
   ],
   "source": [
    "fReader = rawFileReader.fileReader(dir_path+\"data//\"+file_name) # load in the classs object\n",
    "total_chunks = []\n",
    "for start, end, voltage in hv_data.values:\n",
    "    if type(start) == float:\n",
    "        break\n",
    "    print(\"Voltage:\", voltage)\n",
    "    start = \"2024-08-29 \" + start \n",
    "    end = \"2024-08-29 \"+ end\n",
    "    chunks, times, fReader = overview.get_chunks(file_name, fReader=fReader, start=start, end=end)\n",
    "    total_chunks.append((voltage, times, chunks))"
   ]
  },
  {
   "cell_type": "code",
   "execution_count": 23,
   "metadata": {},
   "outputs": [],
   "source": [
    "with open(storage_name, \"wb\") as outp:\n",
    "    pickle.dump(total_chunks, outp)"
   ]
  },
  {
   "cell_type": "code",
   "execution_count": 8,
   "metadata": {},
   "outputs": [
    {
     "name": "stdout",
     "output_type": "stream",
     "text": [
      "Chunks Loaded\n"
     ]
    }
   ],
   "source": [
    "importlib.reload(overview)\n",
    "import matplotlib.pyplot as plt\n",
    "\n",
    "with open(storage_name, \"rb\") as inp:\n",
    "    total_chunks = pickle.load(inp)\n",
    "print(\"Chunks Loaded\")"
   ]
  },
  {
   "cell_type": "code",
   "execution_count": 11,
   "metadata": {},
   "outputs": [
    {
     "name": "stdout",
     "output_type": "stream",
     "text": [
      "# Chunks for 4000: 1\n",
      "Only one chunk, skipping\n",
      "# Chunks for 4200: 1\n",
      "Only one chunk, skipping\n",
      "# Chunks for 4300: 1\n",
      "Only one chunk, skipping\n",
      "# Chunks for 4400: 1\n",
      "Only one chunk, skipping\n",
      "# Chunks for 4500: 1\n",
      "Only one chunk, skipping\n",
      "# Chunks for 4600: 1\n",
      "Only one chunk, skipping\n",
      "# Chunks for 4700: 8\n"
     ]
    },
    {
     "name": "stderr",
     "output_type": "stream",
     "text": [
      "/usera/jd2052/.local/lib/python3.9/site-packages/numpy/_core/fromnumeric.py:3596: RuntimeWarning: Mean of empty slice.\n",
      "  return _methods._mean(a, axis=axis, dtype=dtype,\n",
      "/usera/jd2052/.local/lib/python3.9/site-packages/numpy/_core/_methods.py:140: RuntimeWarning: invalid value encountered in scalar divide\n",
      "  ret = ret / rcount\n"
     ]
    },
    {
     "name": "stdout",
     "output_type": "stream",
     "text": [
      "Cluster Size Done\n",
      "Efficiency Done\n",
      "Possible reconstructions [0, 0, 0, 0, 0, 0]\n",
      "Absolute BVG Done\n",
      "Hit Time Histogram Done\n",
      "# Chunks for 4800: 50\n",
      "Cluster Size Done\n",
      "Efficiency Done\n",
      "Possible reconstructions [0, 0, 0, 0, 0, 0]\n",
      "Absolute BVG Done\n",
      "Hit Time Histogram Done\n",
      "# Chunks for 4900: 164\n",
      "Cluster Size Done\n",
      "Efficiency Done\n",
      "Possible reconstructions [0, 0, 0, 0, 0, 0]\n",
      "Absolute BVG Done\n",
      "Hit Time Histogram Done\n",
      "# Chunks for 5000: 217\n",
      "Cluster Size Done\n",
      "Efficiency Done\n",
      "Possible reconstructions [0, 0, 0, 0, 0, 0]\n",
      "Absolute BVG Done\n",
      "Hit Time Histogram Done\n",
      "# Chunks for 5100: 460\n",
      "Cluster Size Done\n",
      "Efficiency Done\n",
      "Possible reconstructions [0, 0, 0, 0, 0, 0]\n",
      "Absolute BVG Done\n",
      "Hit Time Histogram Done\n",
      "# Chunks for 5200: 455\n",
      "Cluster Size Done\n",
      "Efficiency Done\n",
      "Possible reconstructions [0, 0, 0, 0, 0, 0]\n",
      "Absolute BVG Done\n",
      "Hit Time Histogram Done\n",
      "# Chunks for 5300: 355\n",
      "Cluster Size Done\n",
      "Efficiency Done\n",
      "Possible reconstructions [0, 0, 0, 0, 0, 0]\n",
      "Absolute BVG Done\n",
      "Hit Time Histogram Done\n",
      "# Chunks for 5400: 287\n",
      "Cluster Size Done\n",
      "Efficiency Done\n",
      "Possible reconstructions [0, 0, 0, 0, 0, 0]\n",
      "Absolute BVG Done\n",
      "Hit Time Histogram Done\n",
      "# Chunks for 5450: 418\n",
      "Cluster Size Done\n",
      "Efficiency Done\n",
      "Possible reconstructions [0, 0, 0, 0, 0, 0]\n",
      "Absolute BVG Done\n",
      "Hit Time Histogram Done\n",
      "# Chunks for 5500: 334\n",
      "Cluster Size Done\n",
      "Efficiency Done\n",
      "Possible reconstructions [0, 0, 0, 0, 0, 0]\n",
      "Absolute BVG Done\n",
      "Hit Time Histogram Done\n",
      "# Chunks for 5550: 329\n",
      "Cluster Size Done\n",
      "Efficiency Done\n",
      "Possible reconstructions [0, 0, 0, 0, 0, 0]\n",
      "Absolute BVG Done\n",
      "Hit Time Histogram Done\n",
      "# Chunks for 5600: 267\n",
      "Cluster Size Done\n",
      "Efficiency Done\n",
      "Possible reconstructions [0, 0, 0, 0, 0, 0]\n",
      "Absolute BVG Done\n",
      "Hit Time Histogram Done\n",
      "# Chunks for 5650: 297\n",
      "Cluster Size Done\n",
      "Efficiency Done\n",
      "Possible reconstructions [0, 0, 0, 0, 0, 0]\n",
      "Absolute BVG Done\n",
      "Hit Time Histogram Done\n",
      "# Chunks for 5700: 254\n",
      "Cluster Size Done\n",
      "Efficiency Done\n",
      "Possible reconstructions [0, 0, 0, 0, 0, 0]\n",
      "Absolute BVG Done\n",
      "Hit Time Histogram Done\n",
      "# Chunks for 5750: 359\n",
      "Cluster Size Done\n",
      "Efficiency Done\n",
      "Possible reconstructions [0, 0, 0, 0, 0, 0]\n",
      "Absolute BVG Done\n",
      "Hit Time Histogram Done\n",
      "# Chunks for 5800: 204\n",
      "Cluster Size Done\n",
      "Efficiency Done\n",
      "Possible reconstructions [0, 0, 0, 0, 0, 0]\n",
      "Absolute BVG Done\n",
      "Hit Time Histogram Done\n",
      "# Chunks for 5850: 191\n",
      "Cluster Size Done\n",
      "Efficiency Done\n",
      "Possible reconstructions [0, 0, 0, 0, 0, 0]\n",
      "Absolute BVG Done\n",
      "Hit Time Histogram Done\n",
      "# Chunks for 5900: 194\n",
      "Cluster Size Done\n",
      "Efficiency Done\n",
      "Possible reconstructions [0, 0, 0, 0, 0, 0]\n",
      "Absolute BVG Done\n",
      "Hit Time Histogram Done\n",
      "# Chunks for 5950: 211\n",
      "Cluster Size Done\n",
      "Efficiency Done\n",
      "Possible reconstructions [0, 0, 0, 0, 0, 0]\n",
      "Absolute BVG Done\n",
      "Hit Time Histogram Done\n",
      "# Chunks for 6000: 178\n",
      "Cluster Size Done\n",
      "Efficiency Done\n",
      "Possible reconstructions [0, 0, 0, 0, 0, 0]\n",
      "Absolute BVG Done\n",
      "Hit Time Histogram Done\n"
     ]
    }
   ],
   "source": [
    "importlib.reload(overview)\n",
    "complete_data = []\n",
    "\n",
    "with PdfPages(f\"hit_time_histograms{rpc}.pdf\") as pdf: #trash pdf, so I do not need to click on\n",
    "    for voltage, times, chunks in total_chunks:\n",
    "        results_dict = {}\n",
    "        print(f\"# Chunks for {round(voltage)}:\", len(chunks))\n",
    "        if len(chunks) == 1:\n",
    "            print(\"Only one chunk, skipping\")\n",
    "            continue\n",
    "        \n",
    "        cluster_size, error, hist = overview.cluster_size(chunks, residual=False)\n",
    "        results_dict[\"cluster_size\"] = (cluster_size, error, hist)\n",
    "        results_dict[\"efficiency\"] = overview.efficiency(chunks, residual = False, pdf = pdf)\n",
    "        \n",
    "        \n",
    "        good, bad = overview.abs_bvg_hits(chunks[:100], times[:100], per_rpc=True, pdf = pdf)[rpc]\n",
    "        \n",
    "        results_dict[\"counts\"] = (good, bad)\n",
    "        \n",
    "        hist, std = overview.hit_time_hist(chunks, per_rpc= True, pdf=pdf)\n",
    "        results_dict[\"hit_time_hist\"] = (hist, std)\n",
    "        results_dict[\"voltage\"] = voltage\n",
    "        if results_dict:\n",
    "            complete_data.append(results_dict)\n",
    "\n",
    "\n",
    "        "
   ]
  },
  {
   "cell_type": "code",
   "execution_count": 15,
   "metadata": {},
   "outputs": [
    {
     "name": "stdout",
     "output_type": "stream",
     "text": [
      "{'cluster_size': ([np.float64(nan), np.float64(nan), np.float64(nan), np.float64(nan), np.float64(nan), np.float64(nan)], [100, 100, 100, 100, 100, 100], array([], shape=(6, 0), dtype=object)), 'efficiency': [0, 0, 0, 0, 0, 0], 'counts': ([1740, 365, 106, 451, 804, 769, 1108, 490, 55, 1965, 1531, 1621, 2929, 947, 358, 1282, 804, 72, 643, 312, 847, 1176, 770, 109, 850, 464, 500, 68, 485, 324, 103, 278, 272, 436, 568, 656, 1296, 1605, 4580, 1236, 515, 1173, 273, 460, 919, 2047, 77, 489, 1167, 72, 1451, 1197, 810, 237, 179, 1101, 382, 990, 165, 453, 121, 149, 87, 623, 1366, 315, 154, 590, 328, 564, 722, 974, 382, 118, 264, 199, 419, 396, 1223, 901, 546, 1258, 647, 201, 301, 283, 1419, 676, 412, 331, 372, 709, 261, 1561, 602, 143, 1276, 406, 197, 80], [213, 38, 5, 75, 93, 104, 214, 64, 28, 150, 374, 274, 377, 150, 24, 213, 102, 37, 65, 87, 147, 75, 155, 10, 112, 60, 40, 35, 57, 43, 16, 56, 72, 44, 121, 105, 100, 331, 808, 255, 76, 157, 19, 99, 192, 340, 8, 78, 140, 14, 230, 136, 48, 37, 11, 43, 168, 150, 79, 108, 55, 12, 11, 231, 222, 12, 12, 77, 4, 100, 96, 91, 63, 0, 20, 16, 21, 23, 220, 126, 103, 200, 130, 21, 26, 79, 212, 95, 85, 76, 20, 96, 21, 349, 66, 20, 148, 102, 25, 4]), 'hit_time_hist': ([[(array([   0.  ,   15.99,   31.98,   47.97,   63.96,   79.95,   95.94,\n",
      "        111.93,  127.92,  143.91,  159.9 ,  175.89,  191.88,  207.87,\n",
      "        223.86,  239.85,  255.84,  271.83,  287.82,  303.81,  319.8 ,\n",
      "        335.79,  351.78,  367.77,  383.76,  399.75,  415.74,  431.73,\n",
      "        447.72,  463.71,  479.7 ,  495.69,  511.68,  527.67,  543.66,\n",
      "        559.65,  575.64,  591.63,  607.62,  623.61,  639.6 ,  655.59,\n",
      "        671.58,  687.57,  703.56,  719.55,  735.54,  751.53,  767.52,\n",
      "        783.51,  799.5 ,  815.49,  831.48,  847.47,  863.46,  879.45,\n",
      "        895.44,  911.43,  927.42,  943.41,  959.4 ,  975.39,  991.38,\n",
      "       1007.37, 1023.36, 1039.35, 1055.34, 1071.33, 1087.32, 1103.31,\n",
      "       1119.3 , 1135.29, 1151.28, 1167.27, 1183.26, 1199.25, 1215.24,\n",
      "       1231.23, 1247.22, 1263.21, 1279.2 , 1295.19, 1311.18, 1327.17,\n",
      "       1343.16, 1359.15, 1375.14, 1391.13, 1407.12, 1423.11, 1439.1 ,\n",
      "       1455.09, 1471.08, 1487.07, 1503.06, 1519.05, 1535.04, 1551.03,\n",
      "       1567.02, 1583.01, 1599.  ]), array([1.01133632e-04, 1.07294056e-04, 1.21668379e-04, 1.12427743e-04,\n",
      "       9.49732073e-05, 6.00641365e-05, 2.46416970e-05, 3.08021213e-06,\n",
      "       4.10694950e-06, 9.24063638e-06, 4.41497072e-05, 1.03187106e-04,\n",
      "       1.53497238e-04, 3.61360219e-03, 1.89638393e-02, 2.25096769e-02,\n",
      "       5.89295917e-03, 1.22541106e-03, 7.49004916e-04, 4.47657496e-04,\n",
      "       2.91080046e-04, 1.65818086e-04, 1.63764611e-04, 1.58117556e-04,\n",
      "       1.39636283e-04, 1.51443763e-04, 8.41924648e-05, 1.05240581e-04,\n",
      "       1.07294056e-04, 1.09347531e-04, 1.06780687e-04, 8.36790961e-05,\n",
      "       1.06780687e-04, 1.21155010e-04, 5.69839244e-05, 1.04213844e-04,\n",
      "       8.93261517e-05, 7.75186719e-05, 1.24235222e-04, 8.05988840e-05,\n",
      "       1.08320793e-04, 1.19614904e-04, 9.24063638e-05, 1.05240581e-04,\n",
      "       1.18074798e-04, 8.67593083e-05, 9.39464699e-05, 1.29368909e-04,\n",
      "       1.47850182e-04, 1.37582808e-04, 1.00620263e-04, 7.59785658e-05,\n",
      "       5.18502375e-05, 1.05753950e-04, 8.93261517e-05, 8.16256214e-05,\n",
      "       1.29368909e-04, 9.13796264e-05, 1.36042702e-04, 9.75400507e-05,\n",
      "       8.11122527e-05, 8.57325709e-05, 1.21668379e-04, 6.93047729e-05,\n",
      "       1.16534692e-04, 1.03700475e-04, 1.43229864e-04, 1.09860899e-04,\n",
      "       8.98395204e-05, 8.05988840e-05, 8.98395204e-05, 9.59999446e-05,\n",
      "       1.11914374e-04, 6.87914042e-05, 8.62459396e-05, 5.64705557e-05,\n",
      "       1.06780687e-04, 8.82994143e-05, 9.54865760e-05, 9.54865760e-05,\n",
      "       1.06780687e-04, 1.01647000e-04, 7.23849850e-05, 8.21389901e-05,\n",
      "       6.67379294e-05, 5.28769749e-05, 7.18716163e-05, 9.39464699e-05,\n",
      "       7.18716163e-05, 7.54651971e-05, 1.01647000e-04, 9.29197325e-05,\n",
      "       1.16021323e-04, 1.26802066e-04, 1.15507955e-04, 1.04213844e-04,\n",
      "       9.03528891e-05, 9.13796264e-05, 7.85454092e-05, 6.16042425e-05])), (array([   0.  ,   15.99,   31.98,   47.97,   63.96,   79.95,   95.94,\n",
      "        111.93,  127.92,  143.91,  159.9 ,  175.89,  191.88,  207.87,\n",
      "        223.86,  239.85,  255.84,  271.83,  287.82,  303.81,  319.8 ,\n",
      "        335.79,  351.78,  367.77,  383.76,  399.75,  415.74,  431.73,\n",
      "        447.72,  463.71,  479.7 ,  495.69,  511.68,  527.67,  543.66,\n",
      "        559.65,  575.64,  591.63,  607.62,  623.61,  639.6 ,  655.59,\n",
      "        671.58,  687.57,  703.56,  719.55,  735.54,  751.53,  767.52,\n",
      "        783.51,  799.5 ,  815.49,  831.48,  847.47,  863.46,  879.45,\n",
      "        895.44,  911.43,  927.42,  943.41,  959.4 ,  975.39,  991.38,\n",
      "       1007.37, 1023.36, 1039.35, 1055.34, 1071.33, 1087.32, 1103.31,\n",
      "       1119.3 , 1135.29, 1151.28, 1167.27, 1183.26, 1199.25, 1215.24,\n",
      "       1231.23, 1247.22, 1263.21, 1279.2 , 1295.19, 1311.18, 1327.17,\n",
      "       1343.16, 1359.15, 1375.14, 1391.13, 1407.12, 1423.11, 1439.1 ,\n",
      "       1455.09, 1471.08, 1487.07, 1503.06, 1519.05, 1535.04, 1551.03,\n",
      "       1567.02, 1583.01, 1599.  ]), array([1.23592910e-04, 1.26279712e-04, 1.24667630e-04, 9.40380833e-05,\n",
      "       1.10696258e-04, 5.74975709e-05, 1.88076167e-05, 4.83624428e-06,\n",
      "       9.13512809e-06, 1.28966514e-05, 9.45754438e-05, 1.55834538e-04,\n",
      "       2.18168353e-04, 3.99151362e-03, 2.07055739e-02, 1.99366110e-02,\n",
      "       6.07647226e-03, 1.13114380e-03, 7.79710051e-04, 4.40098230e-04,\n",
      "       2.43961656e-04, 1.45624689e-04, 1.24667630e-04, 1.49923573e-04,\n",
      "       1.30041235e-04, 1.24130270e-04, 1.08009456e-04, 1.07472095e-04,\n",
      "       1.14995142e-04, 1.27891793e-04, 1.39713724e-04, 9.99490486e-05,\n",
      "       9.56501647e-05, 5.53481290e-05, 9.99490486e-05, 1.37026921e-04,\n",
      "       1.12845700e-04, 8.06040714e-05, 7.52304667e-05, 1.04785293e-04,\n",
      "       1.01561130e-04, 1.41863166e-04, 1.05860014e-04, 1.39713724e-04,\n",
      "       1.03710572e-04, 9.24260019e-05, 7.20063038e-05, 9.29633624e-05,\n",
      "       6.77074200e-05, 5.31986871e-05, 8.11414319e-05, 6.44832571e-05,\n",
      "       9.77996066e-05, 1.64432306e-04, 1.19294026e-04, 1.38639003e-04,\n",
      "       8.11414319e-05, 1.26279712e-04, 8.11414319e-05, 9.72622462e-05,\n",
      "       9.67248857e-05, 1.17144584e-04, 1.28966514e-04, 9.45754438e-05,\n",
      "       9.02765600e-05, 8.22161528e-05, 1.24130270e-04, 7.46931062e-05,\n",
      "       9.88743276e-05, 1.18219305e-04, 1.21443468e-04, 5.10492452e-05,\n",
      "       7.95293505e-05, 1.25204991e-04, 1.05860014e-04, 7.95293505e-05,\n",
      "       6.82447805e-05, 8.22161528e-05, 1.10696258e-04, 8.86644786e-05,\n",
      "       1.11770979e-04, 1.29503875e-04, 8.11414319e-05, 8.43655947e-05,\n",
      "       8.22161528e-05, 8.06040714e-05, 9.18886414e-05, 1.26817072e-04,\n",
      "       7.63051876e-05, 6.98568619e-05, 7.30810247e-05, 6.17964548e-05,\n",
      "       8.54403157e-05, 8.38282343e-05, 1.23055549e-04, 1.35952200e-04,\n",
      "       1.20906107e-04, 1.52073015e-04, 1.09084177e-04, 1.08546816e-04])), (array([   0.  ,   15.99,   31.98,   47.97,   63.96,   79.95,   95.94,\n",
      "        111.93,  127.92,  143.91,  159.9 ,  175.89,  191.88,  207.87,\n",
      "        223.86,  239.85,  255.84,  271.83,  287.82,  303.81,  319.8 ,\n",
      "        335.79,  351.78,  367.77,  383.76,  399.75,  415.74,  431.73,\n",
      "        447.72,  463.71,  479.7 ,  495.69,  511.68,  527.67,  543.66,\n",
      "        559.65,  575.64,  591.63,  607.62,  623.61,  639.6 ,  655.59,\n",
      "        671.58,  687.57,  703.56,  719.55,  735.54,  751.53,  767.52,\n",
      "        783.51,  799.5 ,  815.49,  831.48,  847.47,  863.46,  879.45,\n",
      "        895.44,  911.43,  927.42,  943.41,  959.4 ,  975.39,  991.38,\n",
      "       1007.37, 1023.36, 1039.35, 1055.34, 1071.33, 1087.32, 1103.31,\n",
      "       1119.3 , 1135.29, 1151.28, 1167.27, 1183.26, 1199.25, 1215.24,\n",
      "       1231.23, 1247.22, 1263.21, 1279.2 , 1295.19, 1311.18, 1327.17,\n",
      "       1343.16, 1359.15, 1375.14, 1391.13, 1407.12, 1423.11, 1439.1 ,\n",
      "       1455.09, 1471.08, 1487.07, 1503.06, 1519.05, 1535.04, 1551.03,\n",
      "       1567.02, 1583.01, 1599.  ]), array([1.16235569e-04, 1.00074473e-04, 9.57234094e-05, 8.57781201e-05,\n",
      "       9.82097317e-05, 4.16458989e-05, 1.80258368e-05, 2.67279650e-05,\n",
      "       8.08054754e-06, 1.92689980e-05, 8.45349589e-05, 1.83366271e-04,\n",
      "       2.78468100e-04, 4.86013856e-03, 2.21811030e-02, 1.76149721e-02,\n",
      "       5.11436502e-03, 7.35951407e-04, 5.75583617e-04, 3.38761416e-04,\n",
      "       1.88960496e-04, 1.26802438e-04, 1.69069918e-04, 1.36126147e-04,\n",
      "       1.37990889e-04, 9.44802482e-05, 1.04425537e-04, 1.10019763e-04,\n",
      "       1.26802438e-04, 1.26180858e-04, 1.23072955e-04, 1.56016726e-04,\n",
      "       1.46071436e-04, 1.36747728e-04, 1.02560796e-04, 1.19343471e-04,\n",
      "       1.57259887e-04, 1.73420982e-04, 1.50422500e-04, 1.18100310e-04,\n",
      "       9.44802482e-05, 9.88313123e-05, 1.06911860e-04, 1.52908823e-04,\n",
      "       1.52908823e-04, 1.10641343e-04, 1.30531922e-04, 1.23072955e-04,\n",
      "       1.33018244e-04, 1.08155021e-04, 1.07533440e-04, 1.04425537e-04,\n",
      "       1.03182376e-04, 1.20586633e-04, 1.26802438e-04, 1.64718854e-04,\n",
      "       1.34882986e-04, 1.40477211e-04, 1.35504567e-04, 9.07507647e-05,\n",
      "       8.70212812e-05, 1.04425537e-04, 1.16857149e-04, 1.19965052e-04,\n",
      "       1.45449856e-04, 1.12506085e-04, 8.57781201e-05, 8.82644424e-05,\n",
      "       1.18721891e-04, 1.35504567e-04, 1.36747728e-04, 1.31153502e-04,\n",
      "       1.15613988e-04, 1.39234050e-04, 9.88313123e-05, 1.08155021e-04,\n",
      "       8.20486366e-05, 1.27424019e-04, 1.53530403e-04, 9.82097317e-05,\n",
      "       8.76428618e-05, 1.12506085e-04, 1.69069918e-04, 1.51044081e-04,\n",
      "       1.47314598e-04, 1.44828275e-04, 1.19965052e-04, 1.62232531e-04,\n",
      "       1.13127666e-04, 1.42341953e-04, 9.63449900e-05, 1.37990889e-04,\n",
      "       1.33018244e-04, 1.29910341e-04, 1.19343471e-04, 1.12506085e-04,\n",
      "       1.02560796e-04, 1.08776602e-04, 9.63449900e-05, 1.25559277e-04])), (array([   0.  ,   15.99,   31.98,   47.97,   63.96,   79.95,   95.94,\n",
      "        111.93,  127.92,  143.91,  159.9 ,  175.89,  191.88,  207.87,\n",
      "        223.86,  239.85,  255.84,  271.83,  287.82,  303.81,  319.8 ,\n",
      "        335.79,  351.78,  367.77,  383.76,  399.75,  415.74,  431.73,\n",
      "        447.72,  463.71,  479.7 ,  495.69,  511.68,  527.67,  543.66,\n",
      "        559.65,  575.64,  591.63,  607.62,  623.61,  639.6 ,  655.59,\n",
      "        671.58,  687.57,  703.56,  719.55,  735.54,  751.53,  767.52,\n",
      "        783.51,  799.5 ,  815.49,  831.48,  847.47,  863.46,  879.45,\n",
      "        895.44,  911.43,  927.42,  943.41,  959.4 ,  975.39,  991.38,\n",
      "       1007.37, 1023.36, 1039.35, 1055.34, 1071.33, 1087.32, 1103.31,\n",
      "       1119.3 , 1135.29, 1151.28, 1167.27, 1183.26, 1199.25, 1215.24,\n",
      "       1231.23, 1247.22, 1263.21, 1279.2 , 1295.19, 1311.18, 1327.17,\n",
      "       1343.16, 1359.15, 1375.14, 1391.13, 1407.12, 1423.11, 1439.1 ,\n",
      "       1455.09, 1471.08, 1487.07, 1503.06, 1519.05, 1535.04, 1551.03,\n",
      "       1567.02, 1583.01, 1599.  ]), array([7.53996955e-05, 9.97061763e-05, 6.94470880e-05, 5.40695185e-05,\n",
      "       6.69668348e-05, 8.58167587e-05, 1.28973163e-05, 4.46445566e-06,\n",
      "       3.47235440e-06, 1.93459745e-05, 8.73049106e-05, 1.36909973e-04,\n",
      "       1.30957366e-04, 1.53230039e-03, 1.51335126e-02, 2.62728255e-02,\n",
      "       8.64517035e-03, 1.53130829e-03, 9.01820042e-04, 5.86331843e-04,\n",
      "       3.67077465e-04, 1.91475543e-04, 1.33933670e-04, 9.52417207e-05,\n",
      "       9.42496194e-05, 1.07146936e-04, 1.00698278e-04, 9.27614675e-05,\n",
      "       9.57377713e-05, 8.03602018e-05, 8.68088600e-05, 9.97061763e-05,\n",
      "       1.19052151e-04, 6.25023792e-05, 6.64707842e-05, 7.44075943e-05,\n",
      "       7.44075943e-05, 1.06154834e-04, 7.73838980e-05, 6.99431386e-05,\n",
      "       6.00221260e-05, 7.04391892e-05, 7.44075943e-05, 9.17693663e-05,\n",
      "       8.53207081e-05, 7.44075943e-05, 9.77219738e-05, 1.16075847e-04,\n",
      "       9.62338219e-05, 1.00698278e-04, 8.58167587e-05, 8.28404549e-05,\n",
      "       9.62338219e-05, 9.07772650e-05, 6.99431386e-05, 6.99431386e-05,\n",
      "       7.24233918e-05, 1.11611391e-04, 1.19548201e-04, 6.94470880e-05,\n",
      "       6.79589361e-05, 7.14312905e-05, 8.08562524e-05, 7.98641512e-05,\n",
      "       6.89510374e-05, 6.54786830e-05, 5.80379235e-05, 6.05181767e-05,\n",
      "       5.35734679e-05, 7.78799487e-05, 6.74628855e-05, 8.58167587e-05,\n",
      "       8.92891131e-05, 7.93681005e-05, 8.43286068e-05, 8.48246575e-05,\n",
      "       7.53996955e-05, 9.82180244e-05, 9.22654169e-05, 6.25023792e-05,\n",
      "       4.81169110e-05, 6.74628855e-05, 8.18483537e-05, 1.54767796e-04,\n",
      "       1.01690379e-04, 5.80379235e-05, 6.10142273e-05, 1.10619290e-04,\n",
      "       6.74628855e-05, 5.01011135e-05, 9.62338219e-05, 9.57377713e-05,\n",
      "       5.60537210e-05, 7.09352399e-05, 8.82970119e-05, 1.02682480e-04,\n",
      "       1.33437619e-04, 8.53207081e-05, 6.99431386e-05, 5.40695185e-05])), (array([   0.  ,   15.99,   31.98,   47.97,   63.96,   79.95,   95.94,\n",
      "        111.93,  127.92,  143.91,  159.9 ,  175.89,  191.88,  207.87,\n",
      "        223.86,  239.85,  255.84,  271.83,  287.82,  303.81,  319.8 ,\n",
      "        335.79,  351.78,  367.77,  383.76,  399.75,  415.74,  431.73,\n",
      "        447.72,  463.71,  479.7 ,  495.69,  511.68,  527.67,  543.66,\n",
      "        559.65,  575.64,  591.63,  607.62,  623.61,  639.6 ,  655.59,\n",
      "        671.58,  687.57,  703.56,  719.55,  735.54,  751.53,  767.52,\n",
      "        783.51,  799.5 ,  815.49,  831.48,  847.47,  863.46,  879.45,\n",
      "        895.44,  911.43,  927.42,  943.41,  959.4 ,  975.39,  991.38,\n",
      "       1007.37, 1023.36, 1039.35, 1055.34, 1071.33, 1087.32, 1103.31,\n",
      "       1119.3 , 1135.29, 1151.28, 1167.27, 1183.26, 1199.25, 1215.24,\n",
      "       1231.23, 1247.22, 1263.21, 1279.2 , 1295.19, 1311.18, 1327.17,\n",
      "       1343.16, 1359.15, 1375.14, 1391.13, 1407.12, 1423.11, 1439.1 ,\n",
      "       1455.09, 1471.08, 1487.07, 1503.06, 1519.05, 1535.04, 1551.03,\n",
      "       1567.02, 1583.01, 1599.  ]), array([7.34416296e-05, 7.91930825e-05, 9.24656661e-05, 8.18475992e-05,\n",
      "       1.24319867e-04, 1.36265192e-04, 6.37084015e-05, 1.28301642e-05,\n",
      "       1.54846809e-05, 1.41574226e-05, 4.33571066e-05, 1.28744061e-04,\n",
      "       1.18568414e-04, 2.28730858e-04, 5.19798618e-03, 2.23519157e-02,\n",
      "       1.61584858e-02, 4.71486413e-03, 2.75760046e-03, 1.44803888e-03,\n",
      "       7.29549682e-04, 3.04827005e-04, 1.83161654e-04, 1.44228742e-04,\n",
      "       1.06623089e-04, 8.09627603e-05, 8.75990521e-05, 9.91019580e-05,\n",
      "       9.68898607e-05, 7.91930825e-05, 9.64474412e-05, 1.60155843e-04,\n",
      "       1.15913897e-04, 9.37929245e-05, 6.99022739e-05, 7.38840490e-05,\n",
      "       7.43264685e-05, 1.17241156e-04, 1.21222931e-04, 9.82171191e-05,\n",
      "       8.09627603e-05, 8.80414716e-05, 1.11489703e-04, 7.25567907e-05,\n",
      "       1.22550189e-04, 1.46883259e-04, 9.46777634e-05, 8.89263105e-05,\n",
      "       1.18125995e-04, 9.24656661e-05, 5.04358179e-05, 6.63629183e-05,\n",
      "       9.06959883e-05, 1.27416803e-04, 8.62717938e-05, 1.00871636e-04,\n",
      "       1.06180669e-04, 1.38034870e-04, 1.14144219e-04, 1.43343903e-04,\n",
      "       1.22107770e-04, 1.13259381e-04, 1.02641314e-04, 1.01756475e-04,\n",
      "       8.22900187e-05, 1.25204706e-04, 8.27324381e-05, 9.64474412e-05,\n",
      "       7.87506630e-05, 8.00779214e-05, 1.08835186e-04, 1.01756475e-04,\n",
      "       6.23811432e-05, 8.84838910e-05, 8.31748576e-05, 1.00429216e-04,\n",
      "       9.73322801e-05, 1.27416803e-04, 8.75990521e-05, 1.03083733e-04,\n",
      "       9.86595385e-05, 1.00871636e-04, 1.19010833e-04, 1.15913897e-04,\n",
      "       7.07871128e-05, 7.47688879e-05, 6.14963043e-05, 8.62717938e-05,\n",
      "       1.20338092e-04, 1.17241156e-04, 8.18475992e-05, 1.06623089e-04,\n",
      "       1.12374542e-04, 8.00779214e-05, 7.34416296e-05, 5.44175930e-05,\n",
      "       8.75990521e-05, 6.72477572e-05, 9.37929245e-05, 5.48600124e-05])), (array([   0.  ,   15.99,   31.98,   47.97,   63.96,   79.95,   95.94,\n",
      "        111.93,  127.92,  143.91,  159.9 ,  175.89,  191.88,  207.87,\n",
      "        223.86,  239.85,  255.84,  271.83,  287.82,  303.81,  319.8 ,\n",
      "        335.79,  351.78,  367.77,  383.76,  399.75,  415.74,  431.73,\n",
      "        447.72,  463.71,  479.7 ,  495.69,  511.68,  527.67,  543.66,\n",
      "        559.65,  575.64,  591.63,  607.62,  623.61,  639.6 ,  655.59,\n",
      "        671.58,  687.57,  703.56,  719.55,  735.54,  751.53,  767.52,\n",
      "        783.51,  799.5 ,  815.49,  831.48,  847.47,  863.46,  879.45,\n",
      "        895.44,  911.43,  927.42,  943.41,  959.4 ,  975.39,  991.38,\n",
      "       1007.37, 1023.36, 1039.35, 1055.34, 1071.33, 1087.32, 1103.31,\n",
      "       1119.3 , 1135.29, 1151.28, 1167.27, 1183.26, 1199.25, 1215.24,\n",
      "       1231.23, 1247.22, 1263.21, 1279.2 , 1295.19, 1311.18, 1327.17,\n",
      "       1343.16, 1359.15, 1375.14, 1391.13, 1407.12, 1423.11, 1439.1 ,\n",
      "       1455.09, 1471.08, 1487.07, 1503.06, 1519.05, 1535.04, 1551.03,\n",
      "       1567.02, 1583.01, 1599.  ]), array([1.24651532e-04, 1.30220157e-04, 8.52427999e-05, 7.02503477e-05,\n",
      "       7.32488381e-05, 7.66756843e-05, 2.52729909e-05, 8.99547135e-06,\n",
      "       1.02805387e-05, 1.84192985e-05, 7.75323959e-05, 1.44784253e-04,\n",
      "       1.49496167e-04, 8.04023796e-04, 1.37549324e-02, 2.69093100e-02,\n",
      "       8.21243698e-03, 1.30862690e-03, 7.10642236e-04, 5.91987686e-04,\n",
      "       5.62859493e-04, 2.69864140e-04, 1.64060263e-04, 1.44784253e-04,\n",
      "       1.58491638e-04, 1.42642474e-04, 1.13942637e-04, 1.20367974e-04,\n",
      "       1.22509753e-04, 1.28078378e-04, 1.01948675e-04, 1.13514281e-04,\n",
      "       1.41357407e-04, 9.80934732e-05, 1.24651532e-04, 9.42382712e-05,\n",
      "       8.65278672e-05, 1.06660589e-04, 1.10087435e-04, 1.25936599e-04,\n",
      "       1.03662098e-04, 7.71040401e-05, 1.44784253e-04, 1.28935089e-04,\n",
      "       7.58189728e-05, 9.76651175e-05, 1.11372502e-04, 1.17369483e-04,\n",
      "       1.02805387e-04, 1.08374012e-04, 9.08114250e-05, 1.10515791e-04,\n",
      "       1.20796329e-04, 1.41785763e-04, 8.39577326e-05, 9.08114250e-05,\n",
      "       8.43860883e-05, 1.14799349e-04, 1.55064792e-04, 1.06232233e-04,\n",
      "       1.20796329e-04, 1.12657570e-04, 8.73845788e-05, 9.76651175e-05,\n",
      "       1.20796329e-04, 1.14799349e-04, 7.66756843e-05, 1.19511262e-04,\n",
      "       9.93785406e-05, 1.22509753e-04, 1.14370993e-04, 1.32361936e-04,\n",
      "       9.12397808e-05, 1.25508243e-04, 7.36771939e-05, 1.04090454e-04,\n",
      "       1.15227704e-04, 9.59516944e-05, 9.59516944e-05, 1.10087435e-04,\n",
      "       7.79607517e-05, 9.16681366e-05, 9.89501848e-05, 8.35293768e-05,\n",
      "       8.43860883e-05, 9.80934732e-05, 1.11372502e-04, 1.13942637e-04,\n",
      "       1.19082906e-04, 9.46666270e-05, 8.82412903e-05, 9.55233386e-05,\n",
      "       1.18654551e-04, 1.05375521e-04, 1.16512772e-04, 9.38099155e-05,\n",
      "       6.72518572e-05, 7.62473285e-05, 9.38099155e-05, 9.46666270e-05]))]], [[np.float64(265.2985110323465), np.float64(274.36933034879684), np.float64(300.3866963274171), np.float64(249.68915294611412), np.float64(261.36739393233586), np.float64(274.64028876176883)]]), 'voltage': 6000}\n"
     ]
    }
   ],
   "source": [
    "print(complete_data[-1])"
   ]
  },
  {
   "cell_type": "code",
   "execution_count": 12,
   "metadata": {},
   "outputs": [
    {
     "ename": "ValueError",
     "evalue": "max() arg is an empty sequence",
     "output_type": "error",
     "traceback": [
      "\u001b[0;31m---------------------------------------------------------------------------\u001b[0m",
      "\u001b[0;31mValueError\u001b[0m                                Traceback (most recent call last)",
      "Cell \u001b[0;32mIn [12], line 6\u001b[0m\n\u001b[1;32m      4\u001b[0m \u001b[38;5;28;01mfor\u001b[39;00m volt_data \u001b[38;5;129;01min\u001b[39;00m complete_data:\n\u001b[1;32m      5\u001b[0m     \u001b[38;5;28;01mif\u001b[39;00m volt_data[\u001b[38;5;124m\"\u001b[39m\u001b[38;5;124mvoltage\u001b[39m\u001b[38;5;124m\"\u001b[39m] \u001b[38;5;129;01min\u001b[39;00m selected_voltages \u001b[38;5;129;01mand\u001b[39;00m \u001b[38;5;124m\"\u001b[39m\u001b[38;5;124mcluster_size\u001b[39m\u001b[38;5;124m\"\u001b[39m \u001b[38;5;129;01min\u001b[39;00m volt_data\u001b[38;5;241m.\u001b[39mkeys():\n\u001b[0;32m----> 6\u001b[0m         hist, bins \u001b[38;5;241m=\u001b[39m np\u001b[38;5;241m.\u001b[39mhistogram(volt_data[\u001b[38;5;124m\"\u001b[39m\u001b[38;5;124mcluster_size\u001b[39m\u001b[38;5;124m\"\u001b[39m][\u001b[38;5;241m2\u001b[39m][rpc], bins\u001b[38;5;241m=\u001b[39m\u001b[38;5;28mrange\u001b[39m(\u001b[38;5;241m1\u001b[39m, \u001b[38;5;28;43mmax\u001b[39;49m\u001b[43m(\u001b[49m\u001b[43mvolt_data\u001b[49m\u001b[43m[\u001b[49m\u001b[38;5;124;43m\"\u001b[39;49m\u001b[38;5;124;43mcluster_size\u001b[39;49m\u001b[38;5;124;43m\"\u001b[39;49m\u001b[43m]\u001b[49m\u001b[43m[\u001b[49m\u001b[38;5;241;43m2\u001b[39;49m\u001b[43m]\u001b[49m\u001b[43m[\u001b[49m\u001b[43mrpc\u001b[49m\u001b[43m]\u001b[49m\u001b[43m)\u001b[49m\u001b[38;5;241m+\u001b[39m\u001b[38;5;241m2\u001b[39m))\n\u001b[1;32m      7\u001b[0m         plt\u001b[38;5;241m.\u001b[39mstep(bins[:\u001b[38;5;241m-\u001b[39m\u001b[38;5;241m1\u001b[39m], hist, label\u001b[38;5;241m=\u001b[39m\u001b[38;5;124m\"\u001b[39m\u001b[38;5;124mVoltage: \u001b[39m\u001b[38;5;124m\"\u001b[39m\u001b[38;5;241m+\u001b[39mvolt_data[\u001b[38;5;124m\"\u001b[39m\u001b[38;5;124mvoltage\u001b[39m\u001b[38;5;124m\"\u001b[39m]\u001b[38;5;241m+\u001b[39m\u001b[38;5;124m\"\u001b[39m\u001b[38;5;124mV\u001b[39m\u001b[38;5;124m\"\u001b[39m, alpha\u001b[38;5;241m=\u001b[39m\u001b[38;5;241m0.5\u001b[39m)\n\u001b[1;32m      8\u001b[0m plt\u001b[38;5;241m.\u001b[39mlegend()\n",
      "\u001b[0;31mValueError\u001b[0m: max() arg is an empty sequence"
     ]
    }
   ],
   "source": [
    "import matplotlib.pyplot as plt\n",
    "import numpy as np\n",
    "selected_voltages = [5600, 5800, 6000]\n",
    "for volt_data in complete_data:\n",
    "    if volt_data[\"voltage\"] in selected_voltages and \"cluster_size\" in volt_data.keys():\n",
    "        hist, bins = np.histogram(volt_data[\"cluster_size\"][2][rpc], bins=range(1, max(volt_data[\"cluster_size\"][2][rpc])+2))\n",
    "        plt.step(bins[:-1], hist, label=\"Voltage: \"+volt_data[\"voltage\"]+\"V\", alpha=0.5)\n",
    "plt.legend()\n",
    "plt.xlabel(\"Cluster Size\")\n",
    "plt.ylabel(\"Counts\")\n",
    "plt.yscale(\"log\")\n",
    "plt.title(f\"Cluster Size Distribution RPC {rpc}\")\n",
    "plt.show()\n",
    "\n"
   ]
  },
  {
   "cell_type": "code",
   "execution_count": 16,
   "metadata": {},
   "outputs": [
    {
     "name": "stdout",
     "output_type": "stream",
     "text": [
      "Good [3.28, 18.78, 36.45, 45.05, 45.68, 81.54, 120.08, 206.02, 232.52, 228.09, 234.95, 286.25, 257.07, 409.36, 378.03, 422.63, 549.94, 543.93, 541.65, 707.08]\n",
      "Bad [0.41, 1.22, 2.91, 3.45, 3.74, 7.42, 12.18, 23.86, 33.13, 32.67, 31.82, 46.7, 36.0, 64.93, 68.7, 72.77, 103.94, 87.74, 82.86, 107.52]\n",
      "Voltage [4700, 4800, 4900, 5000, 5100, 5200, 5300, 5400, 5450, 5500, 5550, 5600, 5650, 5700, 5750, 5800, 5850, 5900, 5950, 6000]\n"
     ]
    }
   ],
   "source": [
    "import matplotlib.pyplot as plt\n",
    "plt.figure(figsize=(10, 6))  # Create a new figure\n",
    "\n",
    "# Plot In-peak and Off-peak hits in the first figure\n",
    "\n",
    "good = []\n",
    "bad = []\n",
    "voltage = []\n",
    "for volt_data in complete_data:\n",
    "    if \"counts\" not in volt_data.keys():\n",
    "        continue\n",
    "    good.append(sum(volt_data[\"counts\"][0])/100)\n",
    "    bad.append(sum(volt_data[\"counts\"][1])/100)\n",
    "    voltage.append(volt_data[\"voltage\"])\n",
    "\n",
    "plt.plot(voltage, good, label=\"In-peak hits\", color=\"blue\")\n",
    "plt.plot(voltage, bad, label=\"Off-peak hits\", color=\"red\")\n",
    "plt.title(f\"Off-peak and In-peak hits (RPC {rpc})\")\n",
    "plt.xlabel(\"Voltage (V)\")\n",
    "plt.ylabel(\"Counts per chunk\")\n",
    "plt.legend()\n",
    "plt.show()\n",
    "\n",
    "print(\"Good\",good )\n",
    "print(\"Bad\", bad)\n",
    "print(\"Voltage\", voltage)"
   ]
  },
  {
   "cell_type": "code",
   "execution_count": 109,
   "metadata": {},
   "outputs": [],
   "source": [
    "# Plot Efficiency in the third figure\n",
    "plt.figure(figsize=(10, 6))  # Create a new figure\n",
    "efficiencies = []\n",
    "voltage = []\n",
    "for volt_data in complete_data:\n",
    "    if \"efficiency\" not in volt_data.keys():\n",
    "        continue\n",
    "    efficiencies.append(volt_data[\"efficiency\"])\n",
    "    voltage.append(volt_data[\"voltage\"])\n",
    "plt.plot(voltage, efficiencies, label=\"Efficiency\", color=\"red\")\n",
    "plt.title(f\"Efficiency (RPC {rpc})\")\n",
    "plt.xlabel(\"Voltage (V)\")\n",
    "plt.ylabel(\"Efficiency\")\n",
    "plt.legend()\n",
    "plt.show()\n",
    "\n",
    "print(\"Efficiencies:\", efficiencies)\n",
    "print(\"Voltages:\", voltage)"
   ]
  },
  {
   "cell_type": "code",
   "execution_count": 9,
   "metadata": {},
   "outputs": [],
   "source": [
    "# Plot Cluster size in the second figure\n",
    "clusters = []\n",
    "errors = []\n",
    "voltage = []\n",
    "for volt_data in complete_data:\n",
    "    if \"cluster_size\" not in volt_data.keys():\n",
    "        continue\n",
    "    clusters.append(volt_data[\"cluster_size\"][0][rpc])\n",
    "    errors.append(volt_data[\"cluster_size\"][1][rpc])\n",
    "    voltage.append(volt_data[\"voltage\"])\n",
    "\n",
    "plt.figure(figsize=(10, 6))  # Create a new figure\n",
    "plt.plot(voltage, clusters, label=\"Cluster size\", color=\"red\")\n",
    "plt.errorbar(voltage, clusters, yerr=errors, fmt='o', color=\"red\")\n",
    "plt.title(f\"Cluster size (RPC {rpc})\")\n",
    "plt.xlabel(\"Voltage (V)\")\n",
    "plt.ylabel(\"Cluster size\")\n",
    "plt.legend()\n",
    "plt.show()\n",
    "\n",
    "print(\"Cluster size:\", clusters)\n",
    "print(\"Errors:\", errors)\n",
    "print(\"Voltages\", voltage)\n"
   ]
  },
  {
   "cell_type": "code",
   "execution_count": 22,
   "metadata": {},
   "outputs": [],
   "source": [
    "#Plot Histograms\n",
    "rpc = 4\n",
    "selected_voltages = [5600, 5800, 6000]\n",
    "plt.figure(figsize=(10, 6))  # Create a new figure\n",
    "for i, volt_data in enumerate(complete_data):\n",
    "    if volt_data[\"voltage\"] in selected_voltages and \"hit_time_hist\" in volt_data.keys():\n",
    "        hist = volt_data[\"hit_time_hist\"][0][0]\n",
    "        plt.plot(hist[rpc][0][:-1], hist[rpc][1], label=f\"Voltage {volt_data['voltage']} V\")\n",
    "        plt.title(f\"Hit time histogram (RPC {rpc})\")\n",
    "        plt.xlabel(\"Time (25/32 ns)\")\n",
    "        plt.yscale(\"log\")\n",
    "        plt.ylabel(\"Counts (normalised)\")\n",
    "plt.legend()\n",
    "plt.show()"
   ]
  },
  {
   "cell_type": "code",
   "execution_count": 40,
   "metadata": {},
   "outputs": [
    {
     "name": "stdout",
     "output_type": "stream",
     "text": [
      "[[[np.float64(276.95037939735676), np.float64(481.1764271460127), np.float64(261.95351774737145), np.float64(22.580777820499762), np.float64(198.2039152095318), np.float64(105.76022672833302)]], [[np.float64(140.6423798204495), np.float64(194.24480048262524), np.float64(475.3644166274709), np.float64(36.26761181903454), np.float64(240.9350174852225), np.float64(427.37100075479407)]], [[np.float64(422.4884311105226), np.float64(160.47600533105899), np.float64(344.01823110064197), np.float64(115.76890863771862), np.float64(330.8037311703494), np.float64(234.66290423120188)]], [[np.float64(295.6352466000903), np.float64(169.26156145604776), np.float64(286.7519915874869), np.float64(45.15559029300194), np.float64(236.91888095629298), np.float64(274.584223741319)]], [[np.float64(133.89823330118384), np.float64(209.742275299269), np.float64(361.60539299182375), np.float64(15.544845883103225), np.float64(431.5892113587892), np.float64(267.4758099533995)]], [[np.float64(210.66507565369136), np.float64(119.0082447597321), np.float64(319.0290623919765), np.float64(50.23009591777843), np.float64(470.1488982256479), np.float64(363.2432493036059)]], [[np.float64(264.6310217112279), np.float64(270.92739299167374), np.float64(308.6797221105261), np.float64(109.6711771245793), np.float64(197.54606737577475), np.float64(248.54203964127112)]], [[np.float64(265.7226005561734), np.float64(258.3729281231438), np.float64(268.9163488961617), np.float64(117.64960279069595), np.float64(282.190393136297), np.float64(266.027815798506)]], [[np.float64(280.292665701546), np.float64(266.86149834309384), np.float64(275.1214979043538), np.float64(139.9139886519241), np.float64(272.06620787684284), np.float64(266.4196354772713)]], [[np.float64(276.13663203142653), np.float64(264.0487847163466), np.float64(289.63591843500535), np.float64(164.5487685457531), np.float64(270.30126799961886), np.float64(270.0110087595579)]], [[np.float64(266.8334369444558), np.float64(271.039787388424), np.float64(290.34773201311884), np.float64(179.88285706818914), np.float64(266.4560158784681), np.float64(260.67585372011047)]], [[np.float64(271.35944774757417), np.float64(260.9167608182088), np.float64(290.437317274902), np.float64(197.46726684488053), np.float64(272.4049129353581), np.float64(271.1231596090545)]], [[np.float64(269.53263746713833), np.float64(263.7534199825743), np.float64(292.4311780273856), np.float64(204.03595753014383), np.float64(277.15922250620656), np.float64(272.14669876677726)]], [[np.float64(274.12101954826295), np.float64(270.4885779226693), np.float64(288.95772145750016), np.float64(213.83843936997752), np.float64(273.90928027195565), np.float64(274.01928706346285)]], [[np.float64(269.30034237295484), np.float64(267.8741132497974), np.float64(288.40513080913104), np.float64(222.52963942577753), np.float64(276.2660375060531), np.float64(270.33232418491446)]], [[np.float64(271.7642349884049), np.float64(270.55855290549965), np.float64(289.5056630323559), np.float64(228.13897591693993), np.float64(266.33024605203457), np.float64(274.37433492541885)]], [[np.float64(269.61775813356076), np.float64(275.47612929818655), np.float64(293.41318229231445), np.float64(236.74615068110663), np.float64(271.96269712734), np.float64(270.7100470645981)]], [[np.float64(267.83551571213496), np.float64(266.3324950271585), np.float64(288.89723504105916), np.float64(236.3079104812058), np.float64(277.5693787591218), np.float64(275.2814416657471)]], [[np.float64(270.32486439799305), np.float64(267.19716106211223), np.float64(289.4421057884974), np.float64(241.65472667069506), np.float64(267.3947265614691), np.float64(274.1109554223796)]], [[np.float64(273.4434348533184), np.float64(269.3113506855907), np.float64(292.81216157599295), np.float64(252.93947698782003), np.float64(270.19134005652415), np.float64(280.4410334517126)]], [[np.float64(270.08941100310517), np.float64(274.40589734300016), np.float64(294.6241810548917), np.float64(242.48333735286218), np.float64(276.02785532444096), np.float64(270.9922168788531)]], [[np.float64(266.15154753158924), np.float64(272.09716111734133), np.float64(295.6756807006067), np.float64(245.4035123716854), np.float64(273.19587186452566), np.float64(269.16431912968915)]], [[np.float64(269.08397137084944), np.float64(276.6068783639673), np.float64(291.82638512632485), np.float64(241.87816539584296), np.float64(276.0013520655139), np.float64(271.8431363222769)]], [[np.float64(270.74165094635083), np.float64(274.2599982755998), np.float64(287.835551998066), np.float64(240.27277861641534), np.float64(271.1856424187386), np.float64(274.2296139296244)]], [[np.float64(270.2038893979103), np.float64(265.12602863514377), np.float64(287.3202151151581), np.float64(230.55177350620454), np.float64(273.8821132610904), np.float64(269.7460721559368)]], [[np.float64(266.9779156093039), np.float64(270.838693137009), np.float64(289.9119052473732), np.float64(225.56577025801803), np.float64(276.55944801613487), np.float64(271.28479669593236)]]]\n"
     ]
    }
   ],
   "source": [
    "stds = []\n",
    "voltages = []\n",
    "rpc = 3\n",
    "for volt_data in complete_data:\n",
    "    if \"hit_time_hist\" in volt_data.keys():\n",
    "        stds.append(volt_data[\"hit_time_hist\"][1][rpc])\n",
    "        voltages.append(volt_data[\"voltage\"])\n",
    "\n",
    "plt.plot(voltages, stds , label=\"Standard deviation\")\n",
    "plt.xlabel(\"Voltage (V)\")\n",
    "plt.ylabel(\"Standard deviation\")\n",
    "plt.title(f\"Spread of the hit time peak (RPC {rpc})\")\n",
    "plt.legend()\n",
    "plt.show()\n",
    "\n",
    "print(\"Stds:\", stds)\n",
    "print(\"Voltages:\", voltages)"
   ]
  },
  {
   "cell_type": "code",
   "execution_count": 46,
   "metadata": {},
   "outputs": [
    {
     "name": "stdout",
     "output_type": "stream",
     "text": [
      "0\n",
      "1\n",
      "2\n"
     ]
    }
   ],
   "source": [
    "#Save chunks\n",
    "files = [\"proAnubis_240815_1759.raw\", \"proAnubis_240818_1125.raw\", \"proAnubis_240818_1325.raw\"]\n",
    "rpc = 2\n",
    "\n",
    "storage_name = f\"chunks_hv{rpc}.pkl\"\n",
    "file_name = files[rpc]\n",
    "\n",
    "hv_file = 'data/hvScan.csv'  # Replace with your file path\n",
    "hv_data = pd.read_csv(hv_file)\n",
    "#hv_data = pd.read_csv(hv_file, usecols=[\"start_\"+str(rpc),\"end_\"+str(rpc),\"voltage_\"+str(rpc)])\n",
    "\n",
    "#rpc 0: 2024-08-15 17:24:22\t\n",
    "#rpc 1\t2024-08-18 11:39:42\n",
    "#rpc 2\t2024-08-18 12:22:21 \n",
    "\n",
    "\n",
    "import numpy as np\n",
    "import matplotlib.pyplot as plt\n",
    "cl0 = [np.float64(1.0339622641509434), np.float64(1.045432730763078), np.float64(1.0991290983606556), np.float64(1.1098737353933026), np.float64(1.158373989173294), np.float64(1.199706026457619), np.float64(1.243627828510654), np.float64(1.3051052260598126), np.float64(1.3432812085873311), np.float64(1.4434467018370936), np.float64(1.5500040730985418), np.float64(1.6012395661518468), np.float64(1.7499497487437186), np.float64(1.9432747311988978), np.float64(2.1125484678996673), np.float64(2.418690392760291), np.float64(2.816230601735043), np.float64(3.3619210977701544), np.float64(4.264796643667002), np.float64(5.270159270760293)]\n",
    "err0 = [np.float64(0.011126854946426285), np.float64(0.00797484375042157), np.float64(0.010639097472394752), np.float64(0.006399245040924702), np.float64(0.005373077722241347), np.float64(0.011612885000507504), np.float64(0.011525592304870512), np.float64(0.012537685498550555), np.float64(0.011259031916725452), np.float64(0.01603256469828638), np.float64(0.016418504051676275), np.float64(0.016081504890198023), np.float64(0.019236149879349896), np.float64(0.02214773563960034), np.float64(0.026127083422149457), np.float64(0.028637597625313208), np.float64(0.03627862360495955), np.float64(0.03986771546091226), np.float64(0.05496790069215523), np.float64(0.06452445127451842)]\n",
    "eff_0 = [0.7142857142857143, 0.631578947368421, 0.7407407407407407, 0.800498753117207, 0.837734404898584, 0.8513661202185793, 0.8619281045751634, 0.8058705803869246, 0.7712418300653595, 0.7727272727272727, 0.8007478632478633, 0.8371369294605809, 0.06060606060606061, 0.8083756345177665, 0.7622107969151671, 0.6454965357967667, 0.6234747239976758, 0.5983076157292185, 0.5689544579858884, 0.6159267089499648]\n",
    "good_0 = [27.23, 176.75, 187.79, 208.36, 217.25, 219.96, 231.34, 237.86, 244.11, 260.15, 282.66, 307.92, 334.27, 358.79, 387.47, 420.82, 458.28, 558.12, 707.28, 952.21]\n",
    "bad_0 = [0.36, 7.19, 10.49, 13.81, 15.59, 17.85, 20.35, 25.15, 26.74, 27.65, 39.5, 38.18, 45.49, 48.32, 61.55, 60.18, 68.13, 77.08, 101.84, 110.29]\n",
    "std_0 = [np.float64(116.44640183214655), np.float64(147.14398853196343), np.float64(130.4686001108644), np.float64(60.85581179078705), np.float64(147.10844128432203), np.float64(175.56184808696304), np.float64(195.82146789505632), np.float64(203.24078063739753), np.float64(222.7283138879541), np.float64(216.74833196589313), np.float64(233.39542971654942), np.float64(238.1559128032108), np.float64(248.41570803066662), np.float64(254.446597678295), np.float64(254.71945322400353), np.float64(264.04218057542863), np.float64(254.63390691215753), np.float64(260.4522079461875), np.float64(261.20382309648596), np.float64(268.37842991368245), np.float64(265.0101410268362), np.float64(256.88603378265293), np.float64(246.9555944252669)]\n",
    "voltages_0 = [3000.0, 4000.0, 4500.0, 4750.0, 5000.0, 5100.0, 5200.0, 5250.0, 5300.0, 5350.0, 5400.0, 5450.0, 5500.0, 5550.0, 5600.0, 5650.0, 5700.0, 5750.0, 5800.0, 5850.0, 5900.0, 5950.0, 6000.0]\n",
    "\n",
    "\n",
    "cl1 = [np.float64(1.2212389380530972), np.float64(1.1514522821576763), np.float64(1.047808764940239), np.float64(1.074378831881252), np.float64(1.0793083854208945), np.float64(1.130030772943945), np.float64(1.1806598517513687), np.float64(1.2741918045194076), np.float64(1.3286349856293556), np.float64(1.3650258038904328), np.float64(1.4437348498974454), np.float64(1.5276988398916787), np.float64(1.6550848547184367), np.float64(1.7486767031002957), np.float64(1.9628453593196833), np.float64(2.2378623115171603), np.float64(2.5977705099267583), np.float64(3.288475477168398), np.float64(4.1703055084091964), np.float64(5.144818364163138)]\n",
    "err1 = [np.float64(0.15365365868781636), np.float64(0.050306133269046746), np.float64(0.01085851940062489), np.float64(0.012610236466712046), np.float64(0.007588714747234421), np.float64(0.009364968167379703), np.float64(0.010283760210871348), np.float64(0.012014437176205687), np.float64(0.013788835005604861), np.float64(0.012602875350185448), np.float64(0.017094625432900766), np.float64(0.016288766736714727), np.float64(0.01812416969437741), np.float64(0.019781668576846993), np.float64(0.023790335485086862), np.float64(0.029631820011282582), np.float64(0.03495556082164944), np.float64(0.044105284020961893), np.float64(0.05386122345286682), np.float64(0.04904885980725799)]\n",
    "eff_1 = [0, 0.6363636363636364, 0.6, 0.6492890995260664, 0.6837416481069042, 0.7642526964560863, 0.8214285714285714, 0.7619439868204283, 0.7778810408921933, 0.8096013018714402, 0.8154020385050963, 0.8282910874897792, 0.7958030669895076, 0.7861271676300579, 0.7273449920508744, 0.6556451612903226, 0.6175889328063241, 0.6141078838174274, 0.6260543580131209, 0.6921641791044776]\n",
    "good_1 = [23.11, 101.19, 143.96, 161.91, 178.94, 195.78, 209.24, 210.39, 224.92, 246.56, 265.62, 278.61, 294.74, 311.62, 325.12, 350.22, 413.4, 535.85, 680.4, 946.93]\n",
    "bad_1 = [0.82, 2.85, 4.97, 10.8, 11.69, 13.95, 19.5, 22.54, 28.8, 30.71, 39.49, 39.49, 43.2, 48.55, 57.97, 48.03, 68.53, 78.68, 94.05, 143.15]\n",
    "std_1 = [np.float64(33.62848900345012), np.float64(75.40251771826425), np.float64(20.23976986093638), np.float64(98.9403675674452), np.float64(157.65312727971587), np.float64(25.425914261172096), np.float64(155.12154453706253), np.float64(122.58264578206251), np.float64(143.98342398843823), np.float64(176.8360062555137), np.float64(190.7766510634326), np.float64(208.47690585506578), np.float64(231.04675647537343), np.float64(239.5889561940645), np.float64(244.3349516634308), np.float64(242.81730853509387), np.float64(255.74022309532998), np.float64(259.7230028227785), np.float64(263.26304771075047), np.float64(263.73677604056456), np.float64(272.60875767022014), np.float64(264.79901299335796), np.float64(270.115786313314), np.float64(267.52328190045745), np.float64(259.5640546638369), np.float64(256.75948891329324)]\n",
    "voltages_1 = [4000, 4200, 4300, 4400, 4500, 4600, 4700, 4800, 4900, 5000, 5100, 5200, 5300, 5400, 5450, 5500, 5550, 5600, 5650, 5700, 5750, 5800, 5850, 5900, 5950, 6000]\n",
    "\n",
    "cl3 = [np.float64(1.3109243697478992), np.float64(1.4304932735426008), np.float64(1.088797108931337), np.float64(1.0732620320855615), np.float64(1.0889143293264645), np.float64(1.140652965500216), np.float64(1.2034606205250598), np.float64(1.300078919270229), np.float64(1.3564630027632791), np.float64(1.4238995328104238), np.float64(1.585788485880145), np.float64(1.614696409399976), np.float64(1.7720491102649778), np.float64(1.9402696793002916), np.float64(2.200837122032186), np.float64(2.4071046264388887), np.float64(2.8758147771372595), np.float64(3.5219217537402994), np.float64(4.607378273061656), np.float64(5.997637872773203)]\n",
    "err3 = [np.float64(0.13743948558122332), np.float64(0.17042343068733465), np.float64(0.029088649923025792), np.float64(0.014510821401138367), np.float64(0.01037693249897322), np.float64(0.00932813750248753), np.float64(0.01002816662956866), np.float64(0.012962796812449724), np.float64(0.015681357316640058), np.float64(0.015622954630547714), np.float64(0.021914920916389528), np.float64(0.019423017007949815), np.float64(0.02219151199022226), np.float64(0.024554762448446458), np.float64(0.028902655790087), np.float64(0.031205380199562274), np.float64(0.03841785594435168), np.float64(0.047352697168423274), np.float64(0.057293592625362544), np.float64(0.06856915514193454)]\n",
    "eff_3 = [0.3333333333333333, 0, 0.6046511627906976, 0.6592592592592592, 0.7282229965156795, 0.8059701492537313, 0.8470209339774557, 0.8125894134477826, 0.8591352859135286, 0.7796610169491526, 0.7790262172284644, 0.83729662077597, 0.8648960739030023, 0.8652173913043478, 0.8177676537585421, 0.7995712754555199, 0.7423913043478261, 0.6504524886877828, 0.6239035087719298, 0.6424180327868853]\n",
    "good_3 = [39.41, 118.95, 167.0, 178.56, 195.79, 211.84, 229.26, 247.03, 251.18, 267.39, 281.43, 314.74, 349.06, 390.47, 420.62, 477.23, 521.25, 645.4, 786.85, 1058.32]\n",
    "bad_3 =[1.24, 3.28, 5.57, 9.0, 12.72, 14.29, 17.94, 22.83, 22.58, 30.95, 27.29, 30.73, 42.35, 41.41, 49.63, 57.62, 68.54, 78.72, 95.06, 103.75]\n",
    "std_3 = [np.float64(22.580777820499762), np.float64(36.26761181903454), np.float64(115.76890863771862), np.float64(45.15559029300194), np.float64(15.544845883103225), np.float64(50.23009591777843), np.float64(109.6711771245793), np.float64(117.64960279069595), np.float64(139.9139886519241), np.float64(164.5487685457531), np.float64(179.88285706818914), np.float64(197.46726684488053), np.float64(204.03595753014383), np.float64(213.83843936997752), np.float64(222.52963942577753), np.float64(228.13897591693993), np.float64(236.74615068110663), np.float64(236.3079104812058), np.float64(241.65472667069506), np.float64(252.93947698782003), np.float64(242.48333735286218), np.float64(245.4035123716854), np.float64(241.87816539584296), np.float64(240.27277861641534), np.float64(230.55177350620454), np.float64(225.56577025801803)]\n",
    "voltages_3 = [4000, 4200, 4300, 4400, 4500, 4600, 4700, 4800, 4900, 5000, 5100, 5200, 5300, 5400, 5450, 5500, 5550, 5600, 5650, 5700, 5750, 5800, 5850, 5900, 5950, 6000]\n",
    "\n",
    "good_4 = [3.28, 18.78, 36.45, 45.05, 45.68, 81.54, 120.08, 206.02, 232.52, 228.09, 234.95, 286.25, 257.07, 409.36, 378.03, 422.63, 549.94, 543.93, 541.65, 707.08]\n",
    "bad_4 = [0.41, 1.22, 2.91, 3.45, 3.74, 7.42, 12.18, 23.86, 33.13, 32.67, 31.82, 46.7, 36.0, 64.93, 68.7, 72.77, 103.94, 87.74, 82.86, 107.52]\n",
    "voltage_4 = [4700, 4800, 4900, 5000, 5100, 5200, 5300, 5400, 5450, 5500, 5550, 5600, 5650, 5700, 5750, 5800, 5850, 5900, 5950, 6000]\n",
    "\n",
    "\n",
    "\n",
    "cluster = [cl0, cl1, cl3]\n",
    "error = [err0, err1, err3]\n",
    "eff = [eff_0, eff_1, eff_3]\n",
    "good = [good_0, good_1, good_3]\n",
    "bad = [bad_0, bad_1, bad_3]\n",
    "stds = [std_0, std_1, std_3]\n",
    "voltages = [voltages_0, voltages_1, voltages_3]\n",
    "\n",
    "#plot cluster size\n",
    "\"\"\"\n",
    "colour = [\"red\", \"blue\", \"green\"]\n",
    "labels = [\"bottom triplet\", \"middle triplet\", \"singlet\"]\n",
    "plt.figure(figsize=(10, 6))  # Create a new figure\n",
    "for rpc in range(3):\n",
    "    print(rpc)\n",
    "    plt.plot(voltages[rpc][3:], good[rpc], label=labels[rpc]+\" in-peak\", marker = \"o\", color=colour[rpc])\n",
    "    plt.plot(voltages[rpc][3:], bad[rpc], label=labels[rpc]+\" off-peak\", marker = \"x\", linestyle = \"--\", color=colour[rpc])\n",
    "\n",
    "plt.title(\"Hit Counts vs. Voltage\")\n",
    "plt.ylabel(\"hit Count per Chunk\")\n",
    "plt.xlabel(\"Voltage / V\")\n",
    "plt.xlim(5400, 6000)\n",
    "plt.ylim(0, 1200)\n",
    "plt.legend()\n",
    "plt.show()\n",
    "\n",
    "labels = [\"bottom triplet\", \"middle triplet\", \"singlet\"]\n",
    "plt.figure(figsize=(10, 6))  # Create a new figure\n",
    "for rpc in range(3):\n",
    "    print(rpc)\n",
    "    plt.plot(voltages[rpc][3:], cluster[rpc], label=labels[rpc], marker = \"o\", color=colour[rpc])\n",
    "    plt.errorbar(voltages[rpc][3:], cluster[rpc], yerr=error[rpc],  marker = \"o\", color=colour[rpc])\n",
    "plt.title(\"Cluster Size vs. Voltage\")\n",
    "plt.ylabel(\"Average Cluster size\")\n",
    "plt.xlim(5400, 6000)\n",
    "plt.xlabel(\"Voltage / V\")\n",
    "plt.legend()\n",
    "plt.show()\n",
    "\n",
    "labels = [\"bottom triplet\", \"middle triplet\", \"singlet\"]\n",
    "plt.figure(figsize=(10, 6))  # Create a new figure\n",
    "for rpc in range(3):\n",
    "    print(rpc)\n",
    "    plt.plot(voltages[rpc][3:], eff[rpc], label=labels[rpc], marker = \"o\", color=colour[rpc])\n",
    "plt.title(\"Maximal Reconstruction Efficiency vs. Voltage\")\n",
    "plt.ylabel(\"Efficiency\")\n",
    "plt.xlim(5400, 6000)\n",
    "plt.xlabel(\"Voltage / V\")\n",
    "plt.legend()\n",
    "plt.show()\n",
    "\n",
    "\"\"\"\n",
    "\n",
    "labels = [\"bottom triplet\", \"middle triplet\", \"singlet\"]\n",
    "plt.figure(figsize=(10, 6))  # Create a new figure\n",
    "for rpc in range(3):\n",
    "    print(rpc)\n",
    "    plt.plot(voltages[rpc], stds[rpc], label=labels[rpc], marker = \"o\", color=colour[rpc])\n",
    "plt.title(\"Standard deviation of hit time vs. Voltage\")\n",
    "plt.ylabel(\"$\\sigma$ / (25/32) ns\")\n",
    "plt.xlim(5000, 6000)\n",
    "plt.xlabel(\"Voltage / V\")\n",
    "plt.legend()\n",
    "plt.show()\n"
   ]
  }
 ],
 "metadata": {
  "kernelspec": {
   "display_name": "Python 3 (ipykernel)",
   "language": "python",
   "name": "python3"
  },
  "language_info": {
   "codemirror_mode": {
    "name": "ipython",
    "version": 3
   },
   "file_extension": ".py",
   "mimetype": "text/x-python",
   "name": "python",
   "nbconvert_exporter": "python",
   "pygments_lexer": "ipython3",
   "version": "3.9.18"
  }
 },
 "nbformat": 4,
 "nbformat_minor": 4
}
