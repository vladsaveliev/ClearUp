#!/usr/bin/env python
import json
import subprocess
from collections import defaultdict

from geventwebsocket.handler import WebSocketHandler
from gevent.pywsgi import WSGIServer
import gevent

from flask import Flask, jsonify, request, render_template, send_from_directory
from flask_socketio import SocketIO

import time, threading, random, webbrowser, os, platform

from os.path import isfile

from fingerprinting.tree_view import prank_bin

start_local_browser = platform.system() == 'Darwin'

app = Flask(__name__)
socketio = SocketIO(app)

PORT = 5000
MIN_DELAY, MAX_DELAY = 2, 6


time_format = {
    'one': "%H:%M:%S",
    'best': "%a, %d %b %Y %H:%M:%S +0000",
    'other': "%a, %H:%M",
}


# @socketio.on('connected')
# @socketio.on('connected')
# def connected():
#     print('Run_prank, emitting ready')
#     socketio.emit('ready')


@socketio.on('connected')
def run_long_running_tool():
    work_dirpath = '/Users/vlad/vagrant/Fingerprinting/data/Dev_0261_MiSeq_MCRC_PRCC'
    merged_fasta_fpath = '/Users/vlad/vagrant/Fingerprinting/data/Dev_0261_MiSeq_MCRC_PRCC/longprints.best.fas'
    prank_out = os.path.join(work_dirpath, os.path.splitext(os.path.basename(merged_fasta_fpath))[0])

    cmdl = prank_bin + ' -d=' + merged_fasta_fpath + ' -o=' + prank_out + ' -showtree'
    print('Starting prank ' + cmdl)
    proc = subprocess.Popen(cmdl.split(), stderr=subprocess.STDOUT, stdout=subprocess.PIPE)

    lines = []
    prev_time = time.time()
    for stdout_line in iter(proc.stdout.readline, ''):
        print stdout_line.rstrip()
        lines.append(stdout_line)
        cur_time = time.time()
        if cur_time - prev_time > 2:
            socketio.emit('running', json.dumps({
                'finished': False,
                'lines': lines,
            }))
            lines = []
    socketio.emit('running', json.dumps({
        'finished': True,
        'lines': lines,
    }))
    # return render_phylo_tree_page('Dev_0261_MiSeq_MCRC_PRCC')


# @app.route("/data")
# def data():
#     lines = lines_by_runid['']
#     if not lines:
#         raise RuntimeError('Not found lines in ' + str(lines_by_runid))
#
#     print "sending lines " + ''.join(lines_by_runid[''])
#     info = {
#         'lines': lines_by_runid[''],
#         'finished': finished_by_runid[''],
#     }
#     lines_by_runid[''] = []
#     lines_by_runid[''] = None
#     return jsonify(info)


# @app.route("/updated")
# def updated():
#     """
#     Notify the client that an update is ready. Contacted by the client to
#     'subscribe' to the notification service.
#     """
#     ws = request.environ.get('wsgi.websocket', None)
#     print "web socket retrieved"
#     if ws:
#         while True:
#             delay = random.randint(MIN_DELAY, MAX_DELAY)
#             gevent.sleep(delay)
#             ws.send('ready')
#     else:
#         raise RuntimeError("Environment lacks WSGI WebSocket support")


@app.route('/favicon.ico')
def favicon():
    return send_from_directory(os.path.join(app.root_path, 'static'),
                               'favicon.ico', mimetype='image/vnd.microsoft.icon')


@app.route("/")
def main():
    return render_template("serve.html", port=PORT)


if __name__ == "__main__":
    # if start_local_browser:
    #     # start server and web page pointing to it
    #     url = "http://127.0.0.1:{}".format(PORT)
    #     wb = webbrowser.get(None)  # instead of None, can be "firefox" etc
    #     threading.Timer(1.25, lambda: wb.open(url) ).start()
    #
    # http_server = WSGIServer(('', PORT), app, handler_class=WebSocketHandler)
    # http_server.serve_forever()
    # app.run(port=PORT, debug=True)
    socketio.run(app, port=PORT, debug=True)
