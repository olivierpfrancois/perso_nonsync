	def queryEstates(self, layer):
	
		"""
		---------------------------------------------------------------------------------------------
		Function : Opens the query window and handles the query
		---------------------------------------------------------------------------------------------
			Takes 2 arguments:
				-the layer
				-the sql database to query
			The tool allows the user to query the same variable several times, and handles
			the construction of a and or or request.
		---------------------------------------------------------------------------------------------
		"""
		
		#Initialize the Dialog
		query_dialog = QueryDialog()
		
		#Create a dictionary of query variables
		var = {field.name(): field.typeName() for field in layer.pendingFields()}
		for k in var:
			if var[k].lower() == 'varchar':
				var[k] = 'str'
			elif 'int' in var[k].lower():
				var[k] = 'int'
			elif 'float' in var[k].lower():
				var[k] = 'float'
		
		#QMessageBox.warning(self.dlg, "types", '  '.join('{}--{}'.format(key, val) for key, val in var.items()))
		
		#Fill combobox for the area operator
		query_dialog.comboAreaOperator.clear()
		query_dialog.comboAreaOperator.addItems(['=','<','>','<=','>='])
		
		#Get unique values of the Purpose variable
		idx = layer.fieldNameIndex('purpose') #Get index of variable
		purposes = set() #Create empty set for unique values
		for f in layer.getFeatures():
			purposes.add(f.attributes()[idx]) #Populate the set
		
		#Populate the combo boxes for expired and purpose variables
		query_dialog.comboExpired.clear()
		query_dialog.comboExpired.addItems(['All','Yes','No','Unknown'])
		query_dialog.comboPurpose.clear()
		query_dialog.comboPurpose.addItems(['All']+list(purposes))
		
		#Fill comboBox for first extra variable
		query_dialog.comboExtra1Var.clear() #Clear the combo box
		query_dialog.comboExtra1Var.addItems(['-Select-']+sorted(var.keys())) #Add the names to the combo box
		
		#Connect the extra variables combo boxes to filling the next
		query_dialog.comboExtra1Var.currentIndexChanged.connect(lambda: self.populateExtraVar(query_dialog.comboExtra1Var, query_dialog.comboExtra2Var, var))
		query_dialog.comboExtra2Var.currentIndexChanged.connect(lambda: self.populateExtraVar(query_dialog.comboExtra2Var, query_dialog.comboExtra3Var, var))
		query_dialog.comboExtra3Var.currentIndexChanged.connect(lambda: self.populateExtraVar(query_dialog.comboExtra3Var, query_dialog.comboExtra4Var, var))
		
		#Connect the extra variables combo boxes to filling the operator
		query_dialog.comboExtra1Var.currentIndexChanged.connect(lambda: self.populateExtraOp(query_dialog.comboExtra1Var, query_dialog.comboExtra1Operator, var))
		query_dialog.comboExtra2Var.currentIndexChanged.connect(lambda: self.populateExtraOp(query_dialog.comboExtra2Var, query_dialog.comboExtra2Operator, var))
		query_dialog.comboExtra3Var.currentIndexChanged.connect(lambda: self.populateExtraOp(query_dialog.comboExtra3Var, query_dialog.comboExtra3Operator, var))
		query_dialog.comboExtra4Var.currentIndexChanged.connect(lambda: self.populateExtraOp(query_dialog.comboExtra4Var, query_dialog.comboExtra4Operator, var))
				
		#Connect the Cancel button
		query_dialog.cancelButton.clicked.connect(query_dialog.close)
		
		#Connect the apply button
		query_dialog.applyButton.clicked.connect(lambda: self.runQuery(query_dialog, layer, var))
		#Connect the ok button
		query_dialog.okButton.clicked.connect(lambda: self.runQueryOk(query_dialog, layer, var))
		
		#Fix dialog size
		query_dialog.setFixedSize(query_dialog.size())
					
		#Open the Dialog
		query_dialog.show()
		query_dialog.exec_()
	
	
	def populateExtraVar(self, extra_in, extra_out, var):
		if not extra_in.currentText() == '-Select-':
			
			#Clear the combo box
			extra_out.clear()
			
			#Add the names to the combo box
			extra_out.addItems(['-Select-']+sorted(var.keys()))
	
	
	def populateExtraOp(self, extra, extra_op, var):
		#Clear the combo box
		extra_op.clear()
		
		#Get the variable to be used
		v = extra.currentText()
		
		#Add the operators depending on the type of the variable
		if not v == '-Select-':
			if str(var[v]) == 'str':
				extra_op.addItem('Contains')
			else:
				extra_op.addItems(['=','<','>','<=','>='])
	
	
	def getQueryInputs(self, query_dialog, var):
		
		#Create dictionary to hold user inputs
		inputs = {}

		#For each input get the query, the operator and the type
		#The value is a list of lists, with more than one element if the variable was queried twice
		if not query_dialog.lineId.text() == "":
			inputs['idno0'] = [[query_dialog.lineId.text(), 'Contains', var['idno0']]]
		if not query_dialog.lineLessee.text() == "":
			inputs['lesseename'] = [[query_dialog.lineLessee.text(), 'Contains', var['lesseename']]]
		if not query_dialog.lineArea.text() == "":
			inputs['area'] = [[query_dialog.lineArea.text(), query_dialog.comboAreaOperator.currentText(), var['area']]]
		if query_dialog.comboExpired.currentText() == "Yes":
			inputs['expired'] = [["0", '=', var['expired']]]
		elif query_dialog.comboExpired.currentText() == "No":
			inputs['expired'] = [["1", '=', var['expired']]]
		elif query_dialog.comboExpired.currentText() == "Unknown":
			inputs['expired'] = [["2", '=', var['expired']]]
		if not query_dialog.comboPurpose.currentText() == "All":
			inputs['purpose'] = [[query_dialog.comboPurpose.currentText(), 'Contains', var['purpose']]]
		if (not query_dialog.comboExtra1Var.currentText() == "-Select-" and 
				not query_dialog.lineExtra1.text() == ""):
			extra = query_dialog.comboExtra1Var.currentText()
			if not extra in inputs:
				inputs[extra] = [[query_dialog.lineExtra1.text(), 
						query_dialog.comboExtra1Operator.currentText(), 
						var[extra]]]
			else:
				inputs[extra].append([query_dialog.lineExtra1.text(), 
						query_dialog.comboExtra1Operator.currentText(), 
						var[extra]])
			if (not query_dialog.comboExtra2Var.currentText() == "-Select-" and 
					not query_dialog.lineExtra2.text() == ""):
				extra = query_dialog.comboExtra2Var.currentText()
				if not extra in inputs:
					inputs[extra] = [[query_dialog.lineExtra2.text(), 
							query_dialog.comboExtra2Operator.currentText(), 
							var[extra]]]
				else:
					inputs[extra].append([query_dialog.lineExtra2.text(), 
							query_dialog.comboExtra2Operator.currentText(), 
							var[extra]])
				if (not query_dialog.comboExtra3Var.currentText() == "-Select-" and 
						not query_dialog.lineExtra3.text() == ""):
					extra = query_dialog.comboExtra3Var.currentText()
					if not extra in inputs:
						inputs[extra] = [[query_dialog.lineExtra3.text(), 
								query_dialog.comboExtra3Operator.currentText(), 
								var[extra]]]
					else:
						inputs[extra].append([query_dialog.lineExtra3.text(), 
								query_dialog.comboExtra3Operator.currentText(), 
								var[extra]])
					if (not query_dialog.comboExtra4Var.currentText() == "-Select-" and 
							not query_dialog.lineExtra4.text() == ""):
						extra = query_dialog.comboExtra4Var.currentText()
						if not extra in inputs:
							inputs[extra] = [[query_dialog.lineExtra4.text(), 
									query_dialog.comboExtra4Operator.currentText(), 
									var[extra]]]
						else:
							inputs[extra].append([query_dialog.lineExtra4.text(), 
									query_dialog.comboExtra4Operator.currentText(), 
									var[extra]])
		
		#Check if the inputs are correct
		errors = []
		for k, v in inputs.items():
			for cond in v:
				if cond[2] == 'int':
					try:
						int(cond[0])
					except ValueError:
						errors.append(k)
				if cond[2] == 'float':
					try:
						float(cond[0])
					except ValueError:
						errors.append(k)
		
		if len(errors) > 0:
			QMessageBox.critical(self.dlg, "Error Query Inputs", 
					'Values provided for the following variables are wrong'+\
					', '.join(error)+'\n'+\
					'These variables are numerical.')
			return
		else:
			return inputs
		#Gets the inputs from the user and checks that they are valid, or return None
		#Could return a dictionary with single key, 'Error', with the name of the variable(s) for which there is an error
	
	
	def runQuery(self, query_dialog, layer, var):
		#Get user inputs for query
		inputs = self.getQueryInputs(query_dialog, var)
		
		#Exit if inputs are not correct
		if not inputs:
			return
		
		###########################
		#Construct the condition
		Expr = []
		for k,v in inputs.iteritems():
			if v[0][1] == 'Contains':
				expr = ["(\""+k+"\" "+" ILIKE "+" \'%"+cond[0]+"%\'"+" OR "+\
						"\""+k+"\" "+" ILIKE "+" \'"+cond[0]+"%\'"+" OR "+\
						"\""+k+"\" "+" ILIKE "+" \'%"+cond[0]+"\')" for cond in v]
				expr = ' AND '.join(expr)
			
			else:
				#Gather the conditions
				op = [cond[1] for cond in v]
				val = [float(cond[0]) for cond in v]
				e = ["\""+k+"\" "+cond[1]+" "+cond[0] for cond in v]
				
				#Check if there are equal conditions
				equal = [i for i, x in enumerate(op) if x == '=']
				
				#Check if there are less than conditions
				less = [i for i, x in enumerate(op) if x == '<' or x == '<=']
				
				#Check if there are more than conditions
				more = [i for i, x in enumerate(op) if x == '>' or x == '>=']
				
				#Combine the equal conditions using OR
				if len(equal) > 0:
					equal = [e[i] for i in equal]
					equal = '('+' OR '.join(equal)+')'
				
				#If there are only less or more conditions, combine using AND
				if (len(less) > 0 and len(more) == 0):
					diff = [e[i] for i in less]
					diff = '('+' AND '.join(diff)+')'
				elif(len(less) == 0 and len(more) > 0):
					diff = [e[i] for i in more]
					diff = '('+' AND '.join(diff)+')'
				else:
					#If both types are present, keep only the smallest value for less and
					#the biggest value for more
					lessV = [val[i] for i in less]
					lessV = lessV.index(min(lessV))
					lessV = less[lessV]
					moreV = [val[i] for i in more]
					moreV = moreV.index(max(moreV))
					moreV = more[moreV]
					
					if val[lessV] >= val[moreV]:
						#Combine using AND
						diff = '('+e[lessV]+' AND '+e[moreV]+')'
					else:
						#Combine using OR
						diff = '('+e[lessV]+' OR '+e[moreV]+')'
				
				#Combine the equal and less/more conditions using OR
				if len(equal) > 0 and (len(less) > 0 or len(more) > 0):
					expr = '('+equal+' OR '+diff+')'
				elif len(equal) == 0:
					expr = diff
				elif len(less) == 0 and len(more) == 0:
					expr = equal
			
			Expr.append(expr)
		
		#Join the expressions using AND operator
		Expr = ' AND '.join(Expr)
		
		#Check the expression
		#QMessageBox.warning(self.dlg, "expression", Expr)
		
		#Transform into a qgis expression
		Expr = QgsExpression(Expr)
		
		#Get the features
		selec = layer.getFeatures(QgsFeatureRequest(Expr))
		
		#Extract the feature ids
		ids = [i.id() for i in selec]
		
		#Select the features
		layer.setSelectedFeatures(ids)
		
		'''
		#Get the name of the database
		nm = db.databaseName()
		
		#Construct the condition
		sql = "SELECT SYSDEEDID FROM "+nm+" WHERE "
		Expr = []
		for k, v in inputs.items():
			if v[0][1] == 'Contains':
				expr = ["("+k+" LIKE "+" \'%"+cond[0]+"%\'"+" OR "+\
						k+" LIKE "+" \'"+cond[0]+"%\'"+" OR "+\
						k+" LIKE "+" \'%"+cond[0]+"\')" for cond in v]
				expr = ' AND '.join(expr)
			
			else:
				#Gather the conditions
				op = [cond[1] for cond in v]
				val = [float(cond[0]) for cond in v]
				e = [k+" "+cond[1]+" "+cond[0] for cond in v]
				
				#Check if there are equal conditions
				equal = [i for i, x in enumerate(op) if x == '=']
				
				#Check if there are less than conditions
				less = [i for i, x in enumerate(op) if x == '<' or x == '<=']
				
				#Check if there are more than conditions
				more = [i for i, x in enumerate(op) if x == '>' or x == '>=']
				
				#Combine the equal conditions using OR
				if len(equal) > 0:
					equal = [e[i] for i in equal]
					equal = '('+' OR '.join(equal)+')'
				
				#If there are only less or more conditions, combine using AND
				if (len(less) > 0 and len(more) == 0):
					diff = [e[i] for i in less]
					diff = '('+' AND '.join(diff)+')'
				elif(len(less) == 0 and len(more) > 0):
					diff = [e[i] for i in more]
					diff = '('+' AND '.join(diff)+')'
				else:
					#If both types are present, keep only the smallest value for less and
					#the biggest value for more
					lessV = [val[i] for i in less]
					lessV = lessV.index(min(lessV))
					lessV = less[lessV]
					moreV = [val[i] for i in more]
					moreV = moreV.index(max(moreV))
					moreV = more[moreV]
					
					if val[lessV] >= val[moreV]:
						#Combine using AND
						diff = '('+e[lessV]+' AND '+e[moreV]+')'
					else:
						#Combine using OR
						diff = '('+e[lessV]+' OR '+e[moreV]+')'
				
				#Combine the equal and less/more conditions using OR
				if len(equal) > 0 and (len(less) > 0 or len(more) > 0):
					expr = '('+equal+' OR '+diff+')'
				elif len(equal) == 0:
					expr = diff
				elif len(less) == 0 and len(more) == 0:
					expr = equal
			
			Expr.append(expr)
		
		sql = sql + ' AND '.join(Expr) + ';'
		QMessageBox.warning(self.dlg, "expression", sql)
		#Run the query
		query = QSqlQuery(db)
		query.exec_(sql)
		
		#Extract the id into a list
		sql_ids = []
		while query.next(): 
			sql_ids += [query.value(0).toString()]
			#record = query.record()
			#print record.value(0) 
		
		#Get the feature ids for these values in qgis
		'''
	
	def runQueryOk(self, query_dialog, layer, var):
		self.runQuery(query_dialog, layer, var)
		query_dialog.close()
