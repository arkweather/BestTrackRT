import besttrack_RT as btrt
from btengine import btengine
import datetime
import os
import sys
import timeit
import numpy as np
from mpl_toolkits.basemap import Basemap
from multiprocessing import Pool, Manager, Value, Array, Lock
import multiprocessing
from contextlib import closing
import json
import argparse

# Mapping constants
MIN_LAT = 20
MAX_LAT = 51
MIN_LON = -119
MAX_LON = -62

def total_seconds(timedelta):
	return((timedelta.microseconds + 0.0 + (timedelta.seconds + timedelta.days * 24 * 3600) * 10 ** 6) / 10 ** 6)
	
def getOptions():
	parser = argparse.ArgumentParser()
	parser.add_argument('start', type=str, metavar='start', help='Start time in yyyymmdd')
	parser.add_argument('end', type=str, metavar='end', help='End time in yyyymmdd')
	parser.add_argument('inDir', type=str, metavar='inDir', help='Location of source files')
	parser.add_argument('outDir', type=str, metavar='outDir', help='Name of output directory for new tracking files')
	parser.add_argument('-i', '--itype', type=str, metavar='', default='ascii', help='Type of input file: ascii or json')
	parser.add_argument('-t', '--otype', type=str, metavar='', default='ascii', help='Type of output file: ascii, json, or legacy')
	args = parser.parse_args()
	return args

"""
Run the besttrack_RT code for multiple days
"""
def runBTRT(start, end, dates, inDir, outDir, inType, outType):
	
	for thisdate in dates:
	
		dayStormCells = {}
		
		#inPath = '/localdata/ProbSevere/ascii/'
		#inDir = inPath + thisdate.strftime('%Y%m') + '/' + thisdate.strftime('%Y%m%d') + '/'
		times = []
		
		for filename in sorted(os.listdir(inDir)):
			if filename.startswith('._'): continue
			date = str(filename).split('.')[0].split('_')[3]
			time = str(filename).split('.')[0].split('_')[4]
			fileDate = datetime.datetime.strptime(date + '_' + time, '%Y%m%d_%H%M%S')
			print fileDate.date(), thisdate.date()
			if fileDate.date() != thisdate.date(): continue
			times.append(fileDate)
	
		#outPath = '/localdata/ProbSevere/besttrack/btrt/'
		#outDir = outPath + thisdate.strftime('%Y%m') + '/'
		historyPath = outDir + 'history_' + thisdate.strftime('%Y%m%d')+ '.json'
	
		for time in times:
			currentTime = time		
			bufferTime = 3 # minutes
			bufferDist = 10 # km
			historyTime = 30 # minutes 
		
			if outType == 'legacy':	
				stormCells, distanceRatio = btrt.besttrack_RT(currentTime, inDir, '', historyPath, bufferTime, bufferDist, historyTime, outDir, inType, outType)
				dayStormCells.update(stormCells)
			
				bt = btengine(None)
				stormTracks = bt.find_clusters(dayStormCells, dayStormCells.keys())
				stormTracks = bt.theil_sen_batch(stormTracks)
		
				activeCells = dayStormCells.keys()
				generateLegacyOutput(dayStormCells, stormTracks, distanceRatio, outDir, thisdate)
			
			else:
				btrt.besttrack_RT(currentTime, inDir, '', historyPath, bufferTime, bufferDist, historyTime, outDir, inType, outType)
			
		os.remove(historyPath)
		

def generateLegacyOutput(stormCells, stormTracks, distanceRatio, outDir, date):
	
	print 'Preparing output...'
	
	lats = [MIN_LAT, MAX_LAT]
	lons = [MIN_LON, MAX_LON]

	meanLat = np.mean(lats)
	meanLon = np.mean(lons)
	
	m = Basemap(llcrnrlon=MIN_LON, llcrnrlat=MIN_LAT, urcrnrlon=MAX_LON, urcrnrlat=MAX_LAT,
					projection='aeqd', lat_0=meanLat, lon_0=meanLon)
	
	for track in stormTracks:
		
		age = []
		ids = []
		old_IDs = []
	
		for cell in stormTracks[track]['cells']:
			# Convert time back to datetime for processing		
			cell['start_time'] = stormTracks[track]['t0']
			cell['age'] = total_seconds(cell['time'] - cell['start_time'])
			age.append(cell['age'])
			cell['motion_east'] = float(cell['meast'])
			cell['motion_south'] = float(cell['msouth'])
			cell['speed'] = np.sqrt(cell['motion_east'] ** 2 + cell['motion_south'] ** 2)			
			ID = stormCells.keys()[stormCells.values().index(cell)]
			if ID not in ids: ids.append(ID)
			old_IDs.append(cell['oldtrack'])
			
		# Sort cells by age
		ids = [ID for (age, ID) in sorted(zip(age, ids))]			
		
		# Only save cell IDs to storm track to save space
		stormTracks[track]['cells'] = ids
		stormTracks[track]['old_IDs'] = old_IDs

		stormTracks[track]['t0'] = str(stormTracks[track]['t0'])
		stormTracks[track]['tend'] = str(stormTracks[track]['tend'])

		# Convert x, y back to lon, lat and km/s to m/s
		stormTracks[track]['lon0'], stormTracks[track]['lat0'] = m(stormTracks[track]['x0'], stormTracks[track]['y0'], inverse=True)
		stormTracks[track]['lonf'], stormTracks[track]['latf'] = m(stormTracks[track]['xf'], stormTracks[track]['yf'], inverse=True)
		stormTracks[track]['u'] = stormTracks[track]['u'] * distanceRatio * 1000  # m/s
		stormTracks[track]['v'] = stormTracks[track]['v'] * distanceRatio * 1000  # m/s

		# Remove data specific to this run
		stormTracks[track].pop('x0', None)
		stormTracks[track].pop('y0', None)
		stormTracks[track].pop('xf', None)
		stormTracks[track].pop('yf', None)

		if stormTracks.keys().index(track) % 1000 == 0: print '......' + str(stormTracks.keys().index(track)) + ' of ' + str(len(stormTracks)) + ' tracks processed......'
		
	# Finish cleanup
	for cell in stormCells:
		stormCells[cell].pop('x', None)
		stormCells[cell].pop('y', None)
		stormCells[cell].pop('shape_x', None)
		stormCells[cell].pop('shape_y', None)
		stormCells[cell]['time'] = str(stormCells[cell]['time'])
		stormCells[cell]['start_time'] = str(stormCells[cell]['start_time'])
		stormCells[cell].pop('meast', None)
		stormCells[cell].pop('msouth', None)
		
	# Print StormCells to file
	filename = str(date.year) + str(date.month).zfill(2) + str(date.day).zfill(2) + '_cells.data'
	print '\nPrinting ' + filename
	with open(outDir + '/' + filename, 'w') as outfile:
		json.dump(stormCells, outfile, sort_keys=True, indent=0)
		
	# Print StormTracks to file
	filename = str(date.year) + str(date.month).zfill(2) + str(date.day).zfill(2) + '_tracks.data'
	print 'Printing ' + filename
	with open(outDir + '/' + filename, 'w') as outfile:
		json.dump(stormTracks, outfile, sort_keys=True, indent=0)


"""
Main
"""
if __name__ == '__main__':
	
	start_time = timeit.default_timer()

	args = vars(getOptions())
	start = datetime.datetime.strptime(args['start'], '%Y%m%d')
	end = datetime.datetime.strptime(args['end'], '%Y%m%d')
	dates = [start + datetime.timedelta(days=x) for x in range((end-start).days + 1)]
	
	inDir = args['inDir']
	outDir = args['outDir']
	inType = args['itype']
	outType = args['otype']
	
	# TODO: I'm not sure this works...
	if len(dates) > 1:
		subsets = []
		numPerProc = int(np.ceil(float(len(dates)) / multiprocessing.cpu_count()))
		for k in xrange(0, len(dates), numPerProc):
			subsets.append(dates[k:k+numPerProc])
	
		with closing(Pool(processes=multiprocessing.cpu_count(), maxtasksperchild = 1)) as pool:
			[pool.apply_async(runBTRT, (start, end, subsets[l], inDir, outDir, inType, outType,)) for l in range(len(subsets))]
		
			pool.close()
			pool.join()
			pool.terminate()
	
	else:
		runBTRT(start, end, dates, inDir, outDir, inType, outType)

	print timeit.default_timer() - start_time
