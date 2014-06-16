/**********************************
 *** Common routines from k8.js ***
 **********************************/

var getopt = function(args, ostr) {
	var oli; // option letter list index
	if (typeof(getopt.place) == 'undefined')
		getopt.ind = 0, getopt.arg = null, getopt.place = -1;
	if (getopt.place == -1) { // update scanning pointer
		if (getopt.ind >= args.length || args[getopt.ind].charAt(getopt.place = 0) != '-') {
			getopt.place = -1;
			return null;
		}
		if (getopt.place + 1 < args[getopt.ind].length && args[getopt.ind].charAt(++getopt.place) == '-') { // found "--"
			++getopt.ind;
			getopt.place = -1;
			return null;
		}
	}
	var optopt = args[getopt.ind].charAt(getopt.place++); // character checked for validity
	if (optopt == ':' || (oli = ostr.indexOf(optopt)) < 0) {
		if (optopt == '-') return null; //  if the user didn't specify '-' as an option, assume it means null.
		if (getopt.place < 0) ++getopt.ind;
		return '?';
	}
	if (oli+1 >= ostr.length || ostr.charAt(++oli) != ':') { // don't need argument
		getopt.arg = null;
		if (getopt.place < 0 || getopt.place >= args[getopt.ind].length) ++getopt.ind, getopt.place = -1;
	} else { // need an argument
		if (getopt.place >= 0 && getopt.place < args[getopt.ind].length)
			getopt.arg = args[getopt.ind].substr(getopt.place);
		else if (args.length <= ++getopt.ind) { // no arg
			getopt.place = -1;
			if (ostr.length > 0 && ostr.charAt(0) == ':') return ':';
			return '?';
		} else getopt.arg = args[getopt.ind]; // white space
		getopt.place = -1;
		++getopt.ind;
	}
	return optopt;
}

/*************/

function fm_log2tbl(args)
{
	if (args.length == 0) {
		warn("Usage: k8 fermi2.js log2tbl <file1.log> [...]");
		exit(1);
	}
	var buf = new Bytes();
	for (var i = 0; i < args.length; ++i) {
		var fn = args[i];
		var f = new File(fn);
		var prev = 0, curr = 0;
		while (f.readline(buf) >= 0) {
			var m, line = buf.toString();
			if ((m = /mr_restore.*\((\d+)/.exec(line)) != null) {
				prev = m[1];
			} else if ((m = /main_ropebwt2.*\((\d+)/.exec(line)) != null) {
				curr = m[1];
			}
		}
		f.close();
		print(fn, curr, prev);
	}
}

/***********************
 *** Main() function ***
 ***********************/

function main(args)
{
	if (args.length == 0) {
		print("\nUsage:    k8 fermi2.js <command> [arguments]\n");
		print("Commands: log2tbl   extract sample info from ropebwt2 log");
		print("");
		exit(1);
	}

	var cmd = args.shift();
	if (cmd == 'log2tbl') fm_log2tbl(args);
	else warn("Unrecognized command");
}

main(arguments);
