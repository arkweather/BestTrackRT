#from besttrack import btengine
from btengine import btengine
import json
import geojson
import sys
import os
import numpy as np
import datetime
from shapely.geometry.polygon import Polygon
from mpl_toolkits.basemap import Basemap
import scipy.stats.mstats as stats
from bs4 import BeautifulSoup

# Mapping constants
MIN_LAT = 20
MAX_LAT = 51
MIN_LON = -119
MAX_LON = -62

# Function definition
#total_seconds = datetime.timedelta.total_seconds
def total_seconds(timedelta):
	return((timedelta.microseconds + 0.0 + (timedelta.seconds + timedelta.days * 24 * 3600) * 10 ** 6) / 10 ** 6)
	
#==================================================================================================================#
#                                                                                                                  #
#  Read Files                                                                                                      #
#                                                                                                                  #
#==================================================================================================================#

# Handle the old .ascii files from ProbSevere
def readProbSevereAscii(inDir, inSuffix, historyPath, startTime, endTime):
	"""
	Parses probSevere .ascii files
	
	Parameters
	----------
	inDir : string
		The input directory
	inSuffix : string
		The lowest subdirectory of the input directory
	historyPath: string
		The full file path (including extension) of the history json file
	startTime : string
		The earliest time to process
	endTime : string
		The latest time to process
		
	Returns
	-------
	Dictionary
		newCells - Dictionary containing all cells added by the most current file
	Dictionary
		stormCells - Dictionary of storm cells from the history file or past ascii files
	int
		totNumCells - Number of cells in the stormCells dictionary		
	"""
	
	numFiles = 0
	totNumCells = 0
	stormCells = {}
	newCells = {} 
	
	# Try to load the history file
	try:
		print 'Loading storm history...'
		f = open(historyPath)
		stormCells = json.load(f)
		f.close
		
		# Remove cells older than start time
		oldCells = []
		for cell in stormCells:
			stormCells[cell]['time'] = datetime.datetime.strptime(stormCells[cell]['time'], '%Y%m%d_%H%M%S')
			if stormCells[cell]['time'] < startTime: oldCells.append(cell)
		for cell in oldCells:
			stormCells.pop(cell, None)
		
		print 'Number of old cells removed from history: ' + str(len(oldCells))
		
		if len(stormCells) > 0: 
			startTime = endTime
			totNumCells = max([int(key) for key in stormCells.keys()]) + 1
			if totNumCells >= 1e7:
				totNumCells = 0
			print 'Sucessfully loaded history file.  Loading most recent data...'
		else: 
			print 'No recent storms in history.'
			print 'Loading ascii files from ' + str(startTime) + ' to ' + str(endTime) + '...'
			
	# If no history file, make one
	except IOError, ValueError:
		print 'Unable to find storm history file at ' + historyPath + '.'
		print 'Loading ascii files from ' + str(startTime) + ' to ' + str(endTime) + '...'
	
	# Read in ProbSevere files
	for root, dirs, files in os.walk(inDir):
		if inSuffix != '' and not (files and not dirs and os.path.split(root)[-1] == inSuffix): continue
		for asciiFile in files:
			if asciiFile.endswith('.ascii'):
				
				# Skip hidden files
				if asciiFile.startswith('._'): continue
				
				# Check if file falls in date range
				try:
					date = str(asciiFile).split('.')[0].split('_')[3]
					time = str(asciiFile).split('.')[0].split('_')[4]
					fileDate = datetime.datetime.strptime(date + '_' + time, '%Y%m%d_%H%M%S')
				except ValueError:
					print 'File ' + str(asciiFile) + ' has an invalid name.  Expected format SSEC_AWIPS_PROBSEVERE_YYYYMMDD_hhmmss.ascii...'
					continue
				if not startTime <= fileDate <= endTime:
					continue
					
				# Open file
				f = open(root + '/' + asciiFile)
				lines = f.readlines()
				f.close()
				
				print 'Reading ' + asciiFile
				numFiles += 1
				
				
				
				for line in lines:
					if line.startswith('Valid:'): continue
					data = str(line).split(':')
					lats = map(float, data[7].split(',')[0::2])
					lons = map(float, data[7].split(',')[1::2])
					track = data[8]
					prob = int(data[1])
					meast = data[9]
					msouth = data[10]
					
					latr = (max(lats) - min(lats)) / 2.
					lonr = abs(max(lons) - min(lons)) / 2.
					
					# Calculate centroid		
					points = []
					for i in range(0, len(lats)):
						points.append((lons[i], lats[i]))
					poly = Polygon(points)
					
					if not poly.is_valid:
						coords = np.array([float(x) for x in data[7].rstrip().split(',')])
						coords = coords[::-1]
						coords.shape = (len(coords)/2,2)

						polyCoords = [coords.tolist()]
						polyCoords = [[tuple(val) for val in elem] for elem in polyCoords]
						
						polyCoordsCorrected = []
						for coord in polyCoords[0]:
							if coord not in polyCoordsCorrected:
								polyCoordsCorrected.append(coord)
						poly = Polygon(polyCoordsCorrected)

					if not poly.is_valid:
						poly = poly.convex_hull
										
					lon = poly.centroid.x
					lat = poly.centroid.y
					
					cellID = totNumCells
					
					if fileDate == endTime:
						newCells[cellID] = {'prob': prob, 'time': fileDate, 'latr': latr, 'lat': lat, 'lonr': lonr, 'lon': lon, 'meast': meast,
											'msouth': msouth, 'orientation': 'NaN', 'track': track, 'shape_x': lons, 'shape_y': lats, 'ascii': data, 'oldtrack':track}
					else:
						stormCells[cellID] = {'prob': prob, 'time': fileDate, 'latr': latr, 'lat': lat, 'lonr': lonr, 'lon': lon, 'meast': meast,
												'msouth': msouth, 'orientation': 'NaN', 'track': track, 'shape_x': lons, 'shape_y': lats, 'ascii': data, 'oldtrack':track}
					totNumCells += 1
					
	print '\nSuccesfully loaded ' + str(numFiles) + ' ascii files.'
	print 'Number of new cells: ' + str(len(newCells))
	return newCells, stormCells, totNumCells
	
	
# Handle the new json files from ProbSevere	
def readProbSevereJson(inDir, inSuffix, historyPath, startTime, endTime):
	"""
	Parses raw probSevere .json files
	
	Parameters
	----------
	inDir : string
		The input directory
	inSuffix : string
		The lowest subdirectory of the input directory
	historyPath: string
		The full file path (including extension) of the history json file
	startTime : string
		The earliest time to process
	endTime : string
		The latest time to process
		
	Returns
	-------
	Dictionary
		newCells - Dictionary containing all cells added by the most current file
	Dictionary
		stormCells - Dictionary of storm cells from the history file or past ascii files
	int
		totNumCells - Number of cells in the stormCells dictionary		
	"""

	numFiles = 0
	totNumCells = 0
	stormCells = {}
	newCells = {} 
	
	# Try to load the history file
	try:
		print 'Loading storm history...'
		f = open(historyPath)
		stormCells = json.load(f)
		f.close
		
		# Remove cells older than start time
		oldCells = []
		for cell in stormCells:
			stormCells[cell]['time'] = datetime.datetime.strptime(stormCells[cell]['time'], '%Y%m%d_%H%M%S')
			if stormCells[cell]['time'] < startTime: oldCells.append(cell)
		for cell in oldCells:
			stormCells.pop(cell, None)
		
		print 'Number of old cells removed from history: ' + str(len(oldCells))
		
		if len(stormCells) > 0: 
			startTime = endTime
			totNumCells = max([int(key) for key in stormCells.keys()]) + 1
			if totNumCells >= 1e7:
				totNumCells = 0
			print 'Sucessfully loaded history file.  Loading most recent data...'
		else: 
			print 'No recent storms in history.'
			print 'Loading JSON files from ' + str(startTime) + ' to ' + str(endTime) + '...'
			
	# If no history file, make one
	except IOError, ValueError:
		print 'Unable to find storm history file at ' + historyPath + '.'
		print 'Loading ascii files from ' + str(startTime) + ' to ' + str(endTime) + '...'
	
	# Read in ProbSevere files
	for root, dirs, files in os.walk(inDir):
		if inSuffix != '' and not (files and not dirs and os.path.split(root)[-1] == inSuffix): continue
		for jsonFile in files:
			if jsonFile.endswith('.json'):
				
				# Skip hidden files
				if jsonFile.startswith('._'): continue
				
				# Check if file falls in date range
				try:
					date = str(jsonFile).split('.')[0].split('_')[3]
					time = str(jsonFile).split('.')[0].split('_')[4]
					fileDate = datetime.datetime.strptime(date + '_' + time, '%Y%m%d_%H%M%S')
				except ValueError:
					print 'File ' + str(asciiFile) + ' has an invalid name.  Expected format SSEC_AWIPS_PROBSEVERE_YYYYMMDD_hhmmss.json...'
					continue
				if not startTime <= fileDate <= endTime:
					continue
					
				# Open file
				f = open(root + '/' + jsonFile)
				jFile = json.load(f)
				f.close()
				
				print 'Reading ' + jsonFile
				numFiles += 1
				
				index = 0
				for feature in jFile['features']:
					lats = [point[1] for point in feature['geometry']['coordinates'][0]]
					lons = [point[1] for point in feature['geometry']['coordinates'][0]]
					track = int(feature['properties']['ID'])
					prob = int(feature['properties']['PROB'])
					meast = feature['properties']['MOTION_EAST']
					msouth = feature['properties']['MOTION_SOUTH']
					
					latr = (max(lats) - min(lats)) / 2.
					lonr = abs(max(lons) - min(lons)) / 2.
					
					# Calculate centroid
					points = []
					for i in range(0, len(lats)):
						points.append((lons[i], lats[i]))
					poly = Polygon(points)

					# attempt to fix invalid topologies
					if not poly.is_valid:
						polyCoordsCorrected = []
						for coord in feature['geometry']['coordinates'][0]:
							if coord not in polyCoordsCorrected:
								polyCoordsCorrected.append(coord)
						poly = Polygon(polyCoordsCorrected)

					# last ditch attempt *puts face in hands*
					if not poly.is_valid:
						poly = poly.convex_hull
						
					lon = poly.centroid.x
					lat = poly.centroid.y
					
					cellID = totNumCells
					
					if fileDate == endTime:
						newCells[cellID] = {'prob': prob, 'time': fileDate, 'latr': latr, 'lat': lat, 'lonr': lonr, 'lon': lon, 'meast': meast,
											'msouth': msouth, 'orientation': 'NaN', 'track': track, 'shape_x': lons, 'shape_y': lats, 'index': index, 'oldtrack': track}
					else:
						stormCells[cellID] = {'prob': prob, 'time': fileDate, 'latr': latr, 'lat': lat, 'lonr': lonr, 'lon': lon, 'meast':meast,
												'msouth': msouth, 'orientation': 'NaN', 'track': track, 'shape_x': lons, 'shape_y': lats, 'index': index, 'oldtrack':track}
					totNumCells += 1
					index += 1
					
	print '\nSuccesfully loaded ' + str(numFiles) + ' json files.'
	print 'Number of new cells: ' + str(len(newCells))
	return newCells, stormCells, totNumCells, jFile
	
# Handle the xml files from segmotion	
def readSegmotionXML(inDir, inSuffix, historyPath, startTime, endTime):
	"""
	Parses raw segmotion .xml files
	
	Parameters
	----------
	inDir : string
		The input directory
	inSuffix : string
		The lowest subdirectory of the input directory
	historyPath: string
		The full file path (including extension) of the history json file
	startTime : string
		The earliest time to process
	endTime : string
		The latest time to process
		
	Returns
	-------
	Dictionary
		newCells - Dictionary containing all cells added by the most current file
	Dictionary
		stormCells - Dictionary of storm cells from the history file or past ascii files
	int
		totNumCells - Number of cells in the stormCells dictionary		
	"""

	numFiles = 0
	totNumCells = 0
	stormCells = {}
	newCells = {} 
	
	# Try to load the history file
	try:
		print 'Loading storm history...'
		f = open(historyPath)
		stormCells = json.load(f)
		f.close
		
		# Remove cells older than start time
		oldCells = []
		for cell in stormCells:
			stormCells[cell]['time'] = datetime.datetime.strptime(stormCells[cell]['time'], '%Y%m%d_%H%M%S')
			stormCells[cell]['start'] = datetime.datetime.strptime(stormCells[cell]['start'], '%Y%m%d_%H%M%S')
			if stormCells[cell]['time'] < startTime: oldCells.append(cell)
		for cell in oldCells:
			stormCells.pop(cell, None)
		
		print 'Number of old cells removed from history: ' + str(len(oldCells))
		
		if len(stormCells) > 0: 
			startTime = endTime
			totNumCells = max([int(key) for key in stormCells.keys()]) + 1
			if totNumCells >= 1e7:
				totNumCells = 0
			print 'Sucessfully loaded history file.  Loading most recent data...'
		else: 
			print 'No recent storms in history.'
			print 'Loading XML files from ' + str(startTime) + ' to ' + str(endTime) + '...'
			
	# If no history file, make one
	except IOError, ValueError:
		print 'Unable to find storm history file at ' + historyPath + '.'
		print 'Loading XML files from ' + str(startTime) + ' to ' + str(endTime) + '...'
	
	# Read in ProbSevere files
	currentFile = None
	for root, dirs, files in os.walk(inDir):
		if inSuffix != '' and not (files and not dirs and os.path.split(root)[-1] == inSuffix): continue
		for xmlFile in files:
			if xmlFile.endswith('.xml'):
				
				# Skip hidden files
				if xmlFile.startswith('._'): continue
				
				# Check if file falls in date range
				try:
					fileDate = datetime.datetime.strptime(str(xmlFile).split('.')[0], '%Y%m%d-%H%M%S')
				except ValueError:
					print 'File ' + str(xmlFile) + ' has an invalid name.  Expected format YYYYMMDD-hhmmss.xml...'
					continue
				if not startTime <= fileDate <= endTime:
					continue
				if fileDate == endTime: currentFile = xmlFile
				
				# Open file
				f = open(root + '/' + xmlFile)
				xFile = xmlFile
				lines = BeautifulSoup(f, 'html.parser').find_all('datacolumn')
				f.close()
				
				print 'Reading ' + xmlFile
				numFiles += 1
				
				numCells = len(lines[2].find_all('item'))
				
				f = open(root + '/' + xmlFile)
				lines = BeautifulSoup(f, 'html.parser').find_all('data')
				lines = str(lines[0]).split('</datacolumn>')
				f.close()
				
				for i in range(0, numCells):
				
					data = []
					
					for line in lines:			
								
						try:
							name = line.split('"')[1]
							units = line.split('"')[3]
							value = line.split('item value')[i+1].split('"')[1]
							data.append((name, units, value))
						except IndexError:
							continue
				
						if name == 'LatRadius': latr = float(value)
						elif name == 'Latitude': lat = float(value)
						elif name == 'LonRadius': lonr = float(value)
						elif name == 'Longitude': lon = float(value)
						elif name == 'MotionEast': meast = float(value)
						elif name == 'MotionSouth': msouth = float(value)
						elif name == 'Orientation': orientation = float(value)
						elif name == 'RowName': track = int(value)
						elif name == 'Age': age = float(value)
						elif name == 'StartTime': start = datetime.datetime.strptime(value, '%Y%m%d-%H%M%S')
					
					data.append(('OldTrack', 'dimensionless', str(track)))
					cellID = totNumCells
			
					if fileDate == endTime:
						newCells[cellID] = {'prob': 0, 'time': fileDate, 'latr': latr, 'lat': lat, 'lonr': lonr, 'lon': lon, 'meast': meast, 'age': age, 'start': start,
											'msouth': msouth, 'orientation': 'NaN', 'track': track, 'shape_x': [], 'shape_y': [], 'xml': data, 'oldtrack':track}
					else:
						stormCells[cellID] = {'prob': 0, 'time': fileDate, 'latr': latr, 'lat': lat, 'lonr': lonr, 'lon': lon, 'meast': meast, 'age': age, 'start': start,
												'msouth': msouth, 'orientation': orientation, 'track': track, 'shape_x': [], 'shape_y': [], 'xml': data, 'oldtrack':track}
					totNumCells += 1
					
	print '\nSuccesfully loaded ' + str(numFiles) + ' xml files.'
	print 'Number of new cells: ' + str(len(newCells))
	return newCells, stormCells, totNumCells, root + '/' + currentFile

#==================================================================================================================#
#                                                                                                                  #
#  Compare Tracks                                                                                                  #
#                                                                                                                  #
#==================================================================================================================#

def compareTracks(newCells, stormTracks, bufferTime, bufferDist, distanceRatio, existingObjects, modifiedTracksHistory, m, bt):
	"""
	Function to match new cells with existing tracks

	All new cells are compared to each track and 
	added to the most appropriate one based on distance and time

	Parameters
	----------
	newCells : Dictionary
		Dictionary of storCells that don't match an existing track
	stormTracks : Dictionary
		Full stormTracks dictionary containing information about the current 
		tracks and the cells contained within them
	bufferTime : int
		The time threshold to use when associated cells with a track
	bufferDist : int
		The distance threshold to use when associated cells with a track
	distanceRatio : float
		The ratio between x-y distances and lat-lon distances
	existingOjbects : list
		List of new object IDs that are already matched to a track
	modifiedTracksHistory : Dictionary
		Dictionary of previously modified tracks where the key is the original
		object ID and the value is the modified value
	m : Basemap
		Current map projection
	
	Returns
	-------
	Dictionary
		newCells dictionary containing the modified track values
	"""
	
	reservedTracks = []
	changedCells = []
	qlcsObjects = {}
	counter = 0
	qlcsTest = 0

	for cell in newCells:
		# Save cell original track (even if it doesn't change) for visualization
		newCells[cell]['oldtrack'] = newCells[cell]['track']
		
		# Skip new cells that already have a track
		if cell in existingObjects: 
			# Save cell original track (even if it doesn't change) for visualization
			reservedTracks.append(newCells[cell]['track'])
			continue 
		
		cellTime = newCells[cell]['time']
		cellX = newCells[cell]['x']
		cellY = newCells[cell]['y']

		# Calculate distances
		squallRange = []
		minDist = 1e9
		minTrack = newCells[cell]['track']
		
		for track in stormTracks:
			# Make sure two storms don't have the same ID
			if track in reservedTracks: continue			
			# Only compare to tracks in temporal range
			if not (stormTracks[track]['tend'] < cellTime <= stormTracks[track]['tend'] + datetime.timedelta(minutes = bufferTime)):
				continue
			# Avoid having two cells in the same track at the same time
			# The previous if statement should catch this, but explicitly state it here to be sure
			elif stormTracks[track]['tend'] == cellTime: continue

			if stormTracks[track]['u'] == 'NaN':
				xPoint = stormTracks[track]['xf']
				yPoint = stormTracks[track]['yf']
			else:
				xPoint = stormTracks[track]['xf'] + (stormTracks[track]['u'] * (total_seconds(cellTime - stormTracks[track]['tend'])))
				yPoint = stormTracks[track]['yf'] + (stormTracks[track]['v'] * (total_seconds(cellTime - stormTracks[track]['tend'])))

			dist = np.sqrt((cellX - xPoint) ** 2 + (cellY - yPoint) ** 2)
			dist = dist * distanceRatio  # Convert from x,y to km

			if dist < minDist:
				minDist = dist
				minTrack = track
				
			# Get objects in squall line range here to save a step later
			# TODO: Make this range a user setting
			if dist < 150:
				squallRange.append(track)
				
		# If a match is found, replace the original track with the new one and move on
		if minDist <= bufferDist:
			if minTrack != newCells[cell]['track'] and minTrack not in reservedTracks: 
				changedCells.append(cell)
				reservedTracks.append(minTrack)
				newCells[cell]['track'] = minTrack
				stormTracks[minTrack]['cells'].append(newCells[cell])
				stormTracks[minTrack] = bt.theil_sen_single(stormTracks[minTrack])
				
				if counter % 10 == 0:
					print '......' + str(counter) + ' of ' + str(len(newCells) - len(existingObjects)) + ' assigned......'
				counter += 1			
				
				continue
		
		# If no matches found with the conventional algorithm, try the QLCS algorithm
		if len(newCells[cell]['shape_y']) > 0:
			currentlats = map(float, newCells[cell]['shape_y'])
			currentlons = map(float, newCells[cell]['shape_x'])
			currentObj = Polygon([(currentlons[i], currentlats[i]) for i in range(len(currentlats))])
			currentProb = float(newCells[cell]['prob'])
		
			for track in squallRange:
				# Make sure two storms don't have the same ID
				if track in reservedTracks: continue			
			
				# Get last cell in track
				times = {}
				for trackCell in stormTracks[track]['cells']:
					times[trackCell['time'].timetuple()] = trackCell
				lastCell = times[max(times.keys())]
			
				u_vel = (float(lastCell['meast']) / 1000.) / distanceRatio
				v_vel = (-float(lastCell['msouth']) / 1000.) / distanceRatio
			
				# Apply a 5 km buffer around the cell and interpolate to the current time
				oldlats = map(float, lastCell['shape_y'])
				oldlons = map(float, lastCell['shape_x'])
			
				interpPoints = []
			
				# Calculate Centroid
				points = [m(oldlons[i], oldlats[i]) for i in range(len(oldlats))]
				poly = Polygon(points)
			
				if not poly.is_valid:
				
					polyCoords = [points]
					polyCoords = [[tuple(val) for val in elem] for elem in polyCoords]
				
					polyCoordsCorrected = []
					for coord in polyCoords[0]:
						if coord not in polyCoordsCorrected:
							polyCoordsCorrected.append(coord)
					poly = Polygon(polyCoordsCorrected)

				if not poly.is_valid:
					poly = poly.convex_hull

				x0 = poly.centroid.x
				y0 = poly.centroid.y
	
				# Calculate Buffer Points
				for point in points:
					x,y = point
					dx = x - x0
					dy = y - y0
					theta = np.arctan(abs(dy)/float(abs(dx)))
				
					bufferSize = 5. / distanceRatio
					dx2 = abs(bufferSize * np.cos(theta))
					dy2 = abs(bufferSize * np.sin(theta))
		
					if dx < 0: xb = x - dx2
					else: xb = x + dx2
					if dy < 0: yb = y - dy2
					else: yb = y + dy2
				
					# After the buffer is applied, move the object downstream
					xpoint = xb + u_vel * (total_seconds(cellTime - lastCell['time']))
					ypoint = yb + v_vel * (total_seconds(cellTime - lastCell['time']))
		
					interpPoints.append((m(xpoint, ypoint, inverse = True)[0], m(xpoint, ypoint, inverse = True)[1]))
				
				interpObj = Polygon(interpPoints)
			
				# Compare cell to interpolated object
				# Case for splitting QLCS objects
				if interpObj.contains(currentObj.centroid):
					if track not in reservedTracks:
						if track not in qlcsObjects:
							qlcsObjects[track] = {cell: currentProb}
						else:
							qlcsObjects[track][cell] = currentProb
						continue
								
				# Case for merging QLCS objects
				elif currentObj.contains(interpObj.centroid):
					if track not in reservedTracks:
						if track not in qlcsObjects:
							qlcsObjects[track] = {cell: currentProb}
						else:
							qlcsObjects[track][cell] = currentProb
						continue
					
		# Do a final check to see if the object ID has been modified in previous runs
		if newCells[cell]['track'] in modifiedTracksHistory.keys():
			newCells[cell]['track'] = modifiedTracksHistory[newCells[cell]['track']]
			reservedTracks.append(newCells[cell]['track'])
			print '------------History correction----------'
		
		if counter % 10 == 0:
			print '......' + str(counter) + ' of ' + str(len(newCells) - len(existingObjects)) + ' assigned......'
		counter += 1
		
	# Prefer cell with the higher probability for merges and splits
	print 'Handling merges and splits...'
	for track in qlcsObjects:
		maxProb = -1
		preferredCell = -999
		for cell in qlcsObjects[track]:
			if qlcsObjects[track][cell] > maxProb: 
				maxProb = qlcsObjects[track][cell]
				preferredCell = cell
				
		if track != newCells[preferredCell]['track'] and track not in reservedTracks:
			changedCells.append(preferredCell)
			reservedTracks.append(track)
			newCells[preferredCell]['track'] = track
			qlcsTest += 1
			
	print 'QLCS algoritm applied: ' + str(qlcsTest)
	
	return newCells, changedCells
	
#==================================================================================================================#
#                                                                                                                  #
#  Output	                                                                                                       #
#                                                                                                                  #
#==================================================================================================================#	

# Handle old ascii input from ProbSevere		
def outputAscii(currentTime, newCells, stormCells, changedCells, outDir, historyPath):
	"""
	Creates a new ascii and json file with the updated track information
	
	Parameters
	----------
	currentTime : datetime
		The date/time of the most current ascii file
	newCells : Dictionary
		Dictionary of the most current cells
	stormCells : Dictionary
		Dictionary containing all storm cells to be saved
	changedCells : List
		List of track IDs that were changed
	outDir : string
		Filepath where the output files will be saved
	historyPath : string
		The full file path (including extension) of the history json file
	"""
	
	# Save new ascii file
	#filename = 'SSEC_AWIPS_CONVECTPROB_' + currentTime.strftime('%Y%m%d_%H%M%S') + '.ascii'
	filename = 'SSEC_AWIPS_PROBSEVERE_' + currentTime.strftime('%Y%m%d_%H%M%S') + '.ascii'
	print '\nSaving the most recent ascii file: ' + filename
	f = open(outDir + filename, 'w')
	f.write('Valid: ' + currentTime.strftime('%Y%m%d_%H%M%S') + '\n')
	
	for cell in newCells:
		newCells[cell]['ascii'][8] = str(newCells[cell]['track'])
		newCells[cell]['ascii'][-1] = newCells[cell]['ascii'][-1].rstrip()
		newCells[cell]['ascii'].append(str(newCells[cell]['oldtrack']) + '\n')
		#if cell in changedCells: newCells[cell]['ascii'].append('True\n')
		#else: newCells[cell]['ascii'].append('False\n')
		f.write(':'.join(newCells[cell]['ascii']))
	f.close()
	
	# Save new history file
	print 'Saving the new history file...'
	for cell in stormCells:
		stormCells[cell]['time'] = stormCells[cell]['time'].strftime('%Y%m%d_%H%M%S')
	with open(historyPath, 'w') as outfile:
		json.dump(stormCells, outfile, sort_keys = True, indent=1)
	outfile.close()
	
# Handle new json input from ProbSevere	
def outputJson(currentTime, newCells, stormCells, changedCells, outDir, historyPath, jFile):
	"""
	Creates a new json and history json file with the updated track information
	
	Parameters
	----------
	currentTime : datetime
		The date/time of the most current ascii file
	newCells : Dictionary
		Dictionary of the most current cells
	stormCells : Dictionary
		Dictionary containing all storm cells to be saved
	changedCells : List
		List of track IDs that were changed
	outDir : string
		Filepath where the output files will be saved
	historyPath : string
		The full file path (including extension) of the history json file
	"""
	
	# Save new json file
	filename = 'SSEC_AWIPS_PROBSEVERE_' + currentTime.strftime('%Y%m%d_%H%M%S') + '.json'
	print '\nSaving the most recent json file: ' + filename
	
	for cell in newCells:
		index = newCells[cell]['index']
		jFile['features'][index]['properties']['ID'] = newCells[cell]['track']
		jFile['features'][index]['properties']['besttrack'] = newCells[cell]['oldtrack']
		jFile['features'][index]['geometry']['coordinates'] = json.JSONEncoder().encode(jFile['features'][index]['geometry']['coordinates'])

	s = json.dumps(jFile, sort_keys = False, indent = 1).encode('ascii')
	
	with open(outDir + filename, 'w') as outfile:
		for line in s.split('\n'):
			if 'coordinates' in line:
				tmp = line.split(':')
				tmp[1] = tmp[1].replace('"', '')
				line = tmp[0] + ':' + tmp[1]
			outfile.write(line + '\n')
				
	outfile.close()
	
	# Save new history file
	print 'Saving the new history file...'
	for cell in stormCells:
		stormCells[cell]['time'] = stormCells[cell]['time'].strftime('%Y%m%d_%H%M%S')
	with open(historyPath, 'w') as outfile:
		json.dump(stormCells, outfile, sort_keys = True)
	outfile.close()
	
# Handle XML input from segmotion		
def outputXML(currentTime, newCells, stormCells, changedCells, outDir, historyPath, xmlFile):
	"""
	Creates a new xml and json file with the updated track information
	
	Parameters
	----------
	currentTime : datetime
		The date/time of the most current xml file
	newCells : Dictionary
		Dictionary of the most current cells
	stormCells : Dictionary
		Dictionary containing all storm cells to be saved
	changedCells : List
		List of track IDs that were changed
	outDir : string
		Filepath where the output files will be saved
	historyPath : string
		The full file path (including extension) of the history json file
	xmlFile : File
		The current xml file being processed to use as a template
	"""
	
	# Save new xml file
	filename = currentTime.strftime('%Y%m%d-%H%M%S') + '.xml'
	print '\nSaving the most recent xml file: ' + filename
	
	xfile = open(xmlFile)
	lines = xfile.readlines()
	xfile.close()
	
	f = open(outDir + filename, 'w')	
	
	# Don't bother if there are no objects in the file
	if len(newCells.keys()) < 1:
		for line in lines:
			f.write(line)
		f.close()
		
		print 'Saving the new history file...'
		for cell in stormCells:
			stormCells[cell]['time'] = stormCells[cell]['time'].strftime('%Y%m%d_%H%M%S')
		with open(historyPath, 'w') as outfile:
			json.dump(stormCells, outfile, sort_keys = True, indent=1)
		outfile.close()
			
		return
	
	# Otherwise write out the header info
	for line in lines:
		if not line.strip().startswith('<data>'):
			f.write(line)
		else:
			break
	f.write(' <data>\n')
	names = [column[0] for column in newCells[newCells.keys()[0]]['xml']]
	units = [column[1] for column in newCells[newCells.keys()[0]]['xml']]
	values = [[] for i in range(0, len(names))]
	
	for cell in newCells:
		# Last ditch attempt to catch any datetimes that slipped through...
		if type(newCells[cell]['time']) is datetime.datetime:
			newCells[cell]['time'] = newCells[cell]['time'].strftime('%Y%m%d_%H%M%S')
		if type(newCells[cell]['start']) is datetime.datetime:
			newCells[cell]['start'] = newCells[cell]['start'].strftime('%Y%m%d_%H%M%S')
	
		for i in range(0, len(newCells[cell]['xml'])):
			if names[i] == 'RowName': 
				names[i] = 'Track'
				values[i].append(newCells[cell]['track'])
			elif names[i] == 'Age': values[i].append(newCells[cell]['age'])
			elif names[i] == 'StartTime': values[i].append(newCells[cell]['start'])
			elif names[i] == 'MotionEast': values[i].append(str(newCells[cell]['meast']))
			elif names[i] == 'MotionSouth': values[i].append(str(newCells[cell]['msouth']))
			elif names[i] == 'Speed': values[i].append(str(newCells[cell]['speed']))
			else: values[i].append(newCells[cell]['xml'][i][2])
			
	for j in range(0, len(names)):
		f.write('  <datacolumn name="' + names[j] + '" units="' + units[j] + '">\n')
		for k in range(0, len(values[j])):
			f.write('   <item value="' + str(values[j][k]) + '"/>\n')
		f.write('  </datacolumn>\n')
		
	f.write(' </data>\n')
	f.write('</datatable>')
	f.close()
	
	# Save new history file
	print 'Saving the new history file...'
	for cell in stormCells:
		if type(stormCells[cell]['time']) is datetime.datetime:
			stormCells[cell]['time'] = stormCells[cell]['time'].strftime('%Y%m%d_%H%M%S')
		if type(stormCells[cell]['start']) is datetime.datetime:
			stormCells[cell]['start'] = stormCells[cell]['start'].strftime('%Y%m%d_%H%M%S')
	with open(historyPath, 'w') as outfile:
		json.dump(stormCells, outfile, sort_keys = True, indent=1)
	outfile.close()
	
# Convert segmotion xml to json		
def outputSegJson(currentTime, newCells, stormCells, changedCells, outDir, historyPath, xmlFile):
	"""
	Creates a new xml and json file with the updated track information
	
	Parameters
	----------
	currentTime : datetime
		The date/time of the most current xml file
	newCells : Dictionary
		Dictionary of the most current cells
	stormCells : Dictionary
		Dictionary containing all storm cells to be saved
	changedCells : List
		List of track IDs that were changed
	outDir : string
		Filepath where the output files will be saved
	historyPath : string
		The full file path (including extension) of the history json file
	xmlFile : File
		The current xml file being processed to use as a template
	"""
	
	# Save new json file
	dataDicts = []
	
	filename = 'segmotion_' + currentTime.strftime('%Y%m%d-%H%M%S') + '.json'
	print '\nSaving the most recent segmotion-json file: ' + filename
	
	xfile = open(xmlFile)
	lines = xfile.readlines()
	xfile.close()
	
	f = open(outDir + filename, 'w')	
	
	# Don't bother if there are no objects in the file
	if len(newCells.keys()) < 1:
		f.close()
		
		print 'Saving the new history file...'
		for cell in stormCells:
			stormCells[cell]['time'] = stormCells[cell]['time'].strftime('%Y%m%d_%H%M%S')
			stormCells[cell]['start'] = stormCells[cell]['start'].strftime('%Y%m%d_%H%M%S')
		with open(historyPath, 'w') as outfile:
			json.dump(stormCells, outfile, sort_keys = True, indent=1)
		outfile.close()
			
		return
	
	# Otherwise write out the header info
	headerNames = []
	headerValues = []
	
	for line in lines:
		if not line.strip().startswith('<data>'):
			if line.strip().startswith('<datacolumn'):
				headerNames.append(line.strip().split(' ')[1].split('"')[1])
			elif line.strip().startswith('<item'):
				headerValues.append(line.strip().split('"')[1])
		else:
			break
	
	names = [column[0] for column in newCells[newCells.keys()[0]]['xml']]
	
	# Get data from cells
	for cell in newCells:
		
		# Last ditch attempt to catch any datetimes that slipped through...
		if type(newCells[cell]['time']) is datetime.datetime:
			newCells[cell]['time'] = newCells[cell]['time'].strftime('%Y%m%d_%H%M%S')
		if type(newCells[cell]['start']) is datetime.datetime:
			newCells[cell]['start'] = newCells[cell]['start'].strftime('%Y%m%d_%H%M%S')
		
		dataDict = {}
			
		for i in range(len(newCells[cell]['xml'])):
			if names[i] == 'RowName': dataDict['Track'] = int(newCells[cell]['track'])
			elif names[i] == 'Age': dataDict[names[i]] = float(newCells[cell]['age'])
			elif names[i] == 'StartTime': dataDict[names[i]] = newCells[cell]['start']
			elif names[i] == 'OldTrack': dataDict['OldTrack'] = int(newCells[cell]['oldtrack'])
			elif names[i] == 'MotionEast': dataDict['MotionEast'] = float(newCells[cell]['meast'])
			elif names[i] == 'MotionSouth': dataDict['MotionSouth'] = float(newCells[cell]['msouth'])
			elif names[i] == 'Speed': dataDict['Speed'] = float(newCells[cell]['speed'])
			else: dataDict[names[i]] = float(newCells[cell]['xml'][i][2])
			
		for i in range(len(headerNames)): 
			dataDict[headerNames[i]] = headerValues[i]
			
		dataDicts.append(dataDict)
			
	# Save JSON file
	with open(outDir + filename, 'w') as outfile:
		for jdict in dataDicts:
			s = json.dumps(jdict, sort_keys = False).encode('ascii')
			outfile.write(s + '\n')
	outfile.close()		
	
	# Save new history file
	print 'Saving the new history file...'
	for cell in stormCells:
		if type(stormCells[cell]['time']) is datetime.datetime:
			stormCells[cell]['time'] = stormCells[cell]['time'].strftime('%Y%m%d_%H%M%S')
		if type(stormCells[cell]['start']) is datetime.datetime:
			stormCells[cell]['start'] = stormCells[cell]['start'].strftime('%Y%m%d_%H%M%S')
	with open(historyPath, 'w') as outfile:
		json.dump(stormCells, outfile, sort_keys = True, indent=1)
	outfile.close()


	
#==================================================================================================================#
#                                                                                                                  #
#  Main		                                                                                                       #
#                                                                                                                  #
#==================================================================================================================#
		
def besttrack_RT(currentTime, inDir, inSuffix, historyPath, bufferTime, bufferDist, historyTime, outDir, ftype, outtype = ''):
	"""
	Loads current probSevere object and assigns correct track
	
	Parameters
	----------
	currentTime : datetime
		The date/time of the most current ascii file
	inDir : string
		The input directory
	inSuffix : string
		The lowest subdirectory of the input directory
	historyPath : string
		The full file path (including extension) of the history json file
	bufferTime : int
		The time threshold to use when associating cells with a track (minutes)
	bufferDist : int
		The distance threshold to use when associating cells with a track (km)
	historyTime : int
		How long to keep old cells in the history file (minutes)
	outDir : string
		Filepath where the output files will be saved
	ftype : string
		Type of input files to process (ascii or json or xml)
	"""
		
	print 'Running Best Track RT'
		
	if outtype == '': outtype = ftype
	
	# Compute start time
	dt = datetime.timedelta(minutes = historyTime)
	endTime = currentTime
	startTime = currentTime - dt
	
	# Load storm history and most current file
	if ftype == 'ascii': newCells, stormCells, totNumCells = readProbSevereAscii(inDir, inSuffix, historyPath, startTime, endTime)
	elif ftype == 'json': newCells, stormCells, totNumCells, jFile = readProbSevereJson(inDir, inSuffix, historyPath, startTime, endTime)
	elif ftype == 'xml': newCells, stormCells, totNumCells, xmlFile = readSegmotionXML(inDir, inSuffix, historyPath, startTime, endTime)
	else: 
		print 'Invalid file type. Expected "ascii" or "json" or "xml"'
		return
	
	# Set up the btengine
	bt = btengine(None)
	
	#==================================================================================================================#
	#                                                                                                                  #
	#  Map projection                                                                                                  #
	#                                                                                                                  #
	#==================================================================================================================#
	
	# History of changed objects
	modifiedTracksHistory = {}
	
	# Projection variables
	meanLat = np.mean([MIN_LAT, MAX_LAT])
	meanLon = np.mean([MIN_LON, MAX_LON])
	xyDistMax = 0
	llDistMax = 0
	distanceRatio = 0

	# Setup equidistant map projection
	m = Basemap(llcrnrlon=MIN_LON, llcrnrlat=MIN_LAT, urcrnrlon=MAX_LON, urcrnrlat=MAX_LAT,
			    projection='aeqd', lat_0=meanLat, lon_0=meanLon)

	for cell in stormCells:
		stormCells[cell]['x'] = m(stormCells[cell]['lon'], stormCells[cell]['lat'])[0]
		stormCells[cell]['y'] = m(stormCells[cell]['lon'], stormCells[cell]['lat'])[1]
		# Since we're already iterating, get the object's modification history
		modifiedTracksHistory[stormCells[cell]['oldtrack']] = stormCells[cell]['track']
	for cell in newCells:
		newCells[cell]['x'] = m(newCells[cell]['lon'], newCells[cell]['lat'])[0]
		newCells[cell]['y'] = m(newCells[cell]['lon'], newCells[cell]['lat'])[1]

	# Find ratio between x-y distances and lat-lon distances
	xMin, yMin = m(MIN_LON, MIN_LAT)
	xMax, yMax = m(MAX_LON, MAX_LAT)

	xyDistMax = np.sqrt((xMin - xMax) ** 2 + (yMin - yMax) ** 2)

	# Find distance between two lat lon coordinates
	# Source: https://en.wikipedia.org/wiki/Great-circle_distance
	# point1 = [MAX_LON, MIN_LAT]
	# point2 = [MIN_LON, MAX_LAT]

	rlat1 = np.radians(MIN_LAT)
	rlat2 = np.radians(MAX_LAT)
	r = 6371  # Mean radius of Earth (km)
	dlon = abs(MAX_LON - MIN_LON)
	dsig = np.arccos(np.sin(rlat1) * np.sin(rlat2) + np.cos(rlat1) * np.cos(rlat2) * np.cos(np.radians(dlon)))
	llDistMax = r * dsig

	distanceRatio = llDistMax / xyDistMax
	
	#==================================================================================================================#
	#                                                                                                                  #
	#  Cluster Identification                                                                                          #
	#                                                                                                                  #
	#==================================================================================================================#
	
	# Break down into tracks
	print '\nFinding storm history tracks...'
	stormTracks = bt.find_clusters(stormCells, stormCells.keys())
	print 'Number of tracks: ' + str(len(stormTracks))
	
	# Identify potentially bad objects
	print '\nAdding new cells to existing tracks...'
	existingObjects = []
	for cell in newCells:
		# If the cell belongs to an existing track, add it
		if newCells[cell]['track'] in stormTracks:
			stormCells[cell] = newCells[cell]
			existingObjects.append(cell)
	
	print 'Number of cells added to existing tracks: ' + str(len(existingObjects))
	print 'Number of cells with new track ID: ' + str(len(newCells) - len(existingObjects))
	
	#==================================================================================================================#
	#                                                                                                                  #
	#  Track comparison                                                                                                #
	#                                                                                                                  #
	#==================================================================================================================#
	
	# Compare new tracks with existing ones
	print '\nComparing new tracks with existing ones...'
	stormTracks = bt.find_clusters(stormCells, stormCells.keys())
	stormTracks = bt.theil_sen_batch(stormTracks)
	
	newCells, changedCells = compareTracks(newCells, stormTracks, bufferTime, bufferDist, distanceRatio, existingObjects, modifiedTracksHistory, m, bt)
	
	# Update the tracks
	for cell in newCells:
		# If the cell belongs to an existing track, add it
		if newCells[cell]['track'] in stormTracks and cell not in existingObjects:
			stormCells[cell] = newCells[cell]
			existingObjects.append(cell)
	
	print 'Final number of cells added to existing tracks: ' + str(len(existingObjects))		
	print 'Final number of cells with new track ID: ' + str(len(newCells) - len(existingObjects))
	print 'Number of broken tracks fixed: ' + str(len(changedCells))
	
	# Put all cells into the stormCells Dict
	for cell in newCells:
		stormCells[cell] = newCells[cell]
		
	# Update cell age and start times
	# This only works on segmotion data right now
	if ftype == 'xml':
		stormTracks = bt.find_clusters(stormCells, stormCells.keys())
		
		for track in stormTracks:
			times = []
			ids = []
		
			# Sort cells by time
			for cell in stormTracks[track]['cells']:
				ID = stormCells.keys()[stormCells.values().index(cell)]
				if ID not in ids: 
					ids.append(ID)
					times.append(cell['time'])
			
			ids = [ID for (time, ID) in sorted(zip(times, ids))]
		
			for i in range(len(ids)):
				if ids[i] in newCells:
			
					if type(stormCells[ids[0]]['start']) is unicode: # Yes I'm still fighting this stupid issue
						stormCells[ids[0]]['start'] = datetime.datetime.strptime(stormCells[ids[0]]['start'], '%Y%m%d_%H%M%S')
			
					newCells[ids[i]]['start'] = stormCells[ids[0]]['start']
					newCells[ids[i]]['age'] = total_seconds(newCells[ids[i]]['time'] - newCells[ids[i]]['start'])
					if i > 0:
						dt = total_seconds(newCells[ids[i]]['time'] - stormCells[ids[i-1]]['time'])
						dx = newCells[ids[i]]['x'] - stormCells[ids[i-1]]['x']
						dy = newCells[ids[i]]['y'] - stormCells[ids[i-1]]['y']
						try:
							newCells[ids[i]]['meast'] = (dx / dt) * distanceRatio * 1000. # m/s
							newCells[ids[i]]['msouth'] = -(dy / dt) * distanceRatio * 1000. # m/s
						except ZeroDivisionError:
							newCells[ids[i]]['meast'] = 0 # m/s
							newCells[ids[i]]['msouth'] = 0
					newCells[ids[i]]['speed'] = np.sqrt(newCells[ids[i]]['meast']**2 + newCells[ids[i]]['msouth']**2) # m/s

	
	#==================================================================================================================#
	#                                                                                                                  #
	#  Output			                                                                                               #
	#                                                                                                                  #
	#==================================================================================================================#
	
	del stormTracks
	
	# Save output
	if outtype == 'ascii': outputAscii(currentTime, newCells, stormCells, changedCells, outDir, historyPath)
	elif outtype == 'json': outputJson(currentTime, newCells, stormCells, changedCells, outDir, historyPath, jFile)
	elif outtype == 'xml': outputXML(currentTime, newCells, stormCells, changedCells, outDir, historyPath, xmlFile)
	elif outtype == 'seg_json': outputSegJson(currentTime, newCells, stormCells, changedCells, outDir, historyPath, xmlFile)
	
	elif outtype == 'legacy': 
		if not os.path.exists(outDir):
			os.makedirs(outDir)
		print 'Saving the new history file...'
		for cell in stormCells:
			stormCells[cell]['time'] = stormCells[cell]['time'].strftime('%Y%m%d_%H%M%S')
			#stormCells[cell]['start'] = stormCells[ids[i]]['start'].strftime('%Y%m%d_%H%M%S')
		with open(historyPath, 'w') as outfile:
			json.dump(stormCells, outfile, sort_keys = True)
		outfile.close()
		for cell in stormCells:
			stormCells[cell]['time'] = datetime.datetime.strptime(stormCells[cell]['time'],'%Y%m%d_%H%M%S')
			#stormCells[cell]['start'] = datetime.datetime.strptime(stormCells[cell]['start'],'%Y%m%d_%H%M%S')
			
		return stormCells, distanceRatio
	
	else: print 'Something went horribly wrong... (Bad file type)'
	
	print 'Best Track RT complete!\n'
	
	
# Usage example and testing
if __name__ == '__main__':
	
	currentTime = datetime.datetime.strptime('20150506_203639', '%Y%m%d_%H%M%S')
	inDir = '/localdata/ProbSevere/new_json/20160901'
	outDir = '/localdata/ProbSevere/new_json/test2/'
	historyPath = '/localdata/ProbSevere/new_json/test2/history.json'
	bufferTime = 3 # minutes
	bufferDist = 10 # km
	historyTime = 2 * 60 # minutes (2 hours)
	
	besttrack_RT(currentTime, inDir, '', historyPath, bufferTime, bufferDist, historyTime, outDir, 'json')
	
	
	

