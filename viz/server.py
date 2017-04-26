#!/usr/bin/env python

# Load required modules
import sys, os, json, argparse
import tornado.web, tornado.ioloop
from tornado.escape import json_encode

################################################################################
# Index page
class MainHandler(tornado.web.RequestHandler):
	def initialize(self, results):
		self.results = results

	def get(self):
		self.render("index.html", results=self.results)

# Run page
class RunHandler(tornado.web.RequestHandler):
	def initialize(self, result):
		self.run_index = result['run_index']
		self.is_consensus = result['is_consensus']
		self.network_name = result['network_name']
		self.heat_name = result['heat_name']

	def get(self):
		self.render("subnetworks.html", run_index=self.run_index, is_consensus=self.is_consensus, network_name=self.network_name, heat_name=self.heat_name)

# Data pages
class DataHandler(tornado.web.RequestHandler):
	def initialize(self, data):
		self.data = data

	def get(self):
		self.write(json.dumps(self.data))

################################################################################
# Argument parser
def get_parser():
	parser = argparse.ArgumentParser()
	parser.add_argument('-i', '--input_directory', type=str, required=True)
	parser.add_argument('-p', '--port', type=int, required=False, default=8000)
	return parser

def run( args ):
	# Load the results and create a summary for each file
	results = []
	summary = []
	for root in next(os.walk(args.input_directory))[1]:
		# Every subdirectory should have a viz-data.json file, otherwise
		# we skip it
		try:
			with open('{}/{}/viz-data.json'.format(args.input_directory, root), 'r') as IN:
				# Parse the data and extract the info we need
				data = json.load(IN)
				params = data['params']
				is_consensus = params['consensus']
				if is_consensus:
					heat_name, network_name = '', ''
				else:
					heat_name, network_name = params['heat_name'], params['network_name']

				num_subnetworks = len(data['subnetworks'][params['auto_delta']])
				result = dict(heat_name=heat_name, network_name=network_name,
							  auto_delta=params['auto_delta'],
							  num_subnetworks=num_subnetworks, data=data,
							  is_consensus=is_consensus)

				results.append( result )
		except IOError:
			continue

	# Sort so the consensus is first
	results.sort(key=lambda r: not r['is_consensus'])

	# Set up the server
	routes = [ (r'/bower_components/(.*)', tornado.web.StaticFileHandler, {'path': "bower_components"})]
	for i, result in enumerate(results):
		data_route = r"/data/{}.json".format(i)
		routes.append( (data_route, DataHandler, dict(data=result['data'])) )
		run_route = r"/{}".format(i)
		result['run_index'] = i
		routes.append( (run_route, RunHandler, dict(result=result)) )
	routes.append( (r"/", MainHandler, dict(results=results)) )

	# Start server
	app = tornado.web.Application(routes)
	app.listen(args.port)
	print 'Listening on port {}'.format(args.port)
	tornado.ioloop.IOLoop.current().start()

if __name__ == '__main__': run( get_parser().parse_args(sys.argv[1:]) )
