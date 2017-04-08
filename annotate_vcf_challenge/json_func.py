def expand_json_var(json_data, var: str):
	# string of format parent.child.grandchild split into keys for json dict
	levels = var.split(".")
	curr = json_data
	# index into successive dicts
	for lvl in levels:
		if lvl in curr:
			curr = curr[lvl]
		else:
			# if target key doesn't exist, return None, caller will handle condition
			return None
	# after all levels are traversed, we have the value we need
	return curr