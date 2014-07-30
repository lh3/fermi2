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

function fm_id2sam(args)
{
	if (args.length == 0) {
		warn("Usage: k8 fermi2.js id2sam <sample.tbl> <matches.txt>");
		exit(1);
	}
	var buf = new Bytes();
	var f = new File(args[0]);
	var sample = [];
	while (f.readline(buf) >= 0) {
		var t = buf.toString().split("\t");
		sample.push([parseInt(t[1]), t[0]]);
	}
	f.close();
	sample.sort(function(a,b) {return a[0]-b[0];})
	f = args.length > 1? new File(args[1]) : new File();
	while (f.readline(buf) >= 0) {
		var t = buf.toString().split("\t");
		if (t[0] == 'EM') {
			var ids = [];
			for (var i = 4; i < t.length; ++i) {
				var s = t[i].split(":");
				var id = parseInt(s[0]);
				if (id >= sample[sample.length-1][0])
					throw Error("Invalid sequence ID "+id);
				ids.push(id);
			}
			ids.sort(function(a,b){return a-b;});
			for (var i = 0; i < ids.length; ++i) {
				var id = ids[i];
				var start = 0, end = sample.length;
				while (start < end) {
					var mid = (start + end) >> 1;
					if (sample[mid][0] == id) start = end = mid - 1;
					else if (sample[mid][0] < id) start = mid + 1;
					else end = mid;
				}
				t[i+4] = sample[start][1];
			}
			print(t.join("\t"));
		} else print(buf);
	}
	f.close();
}

function fm_cnt2hist(args)
{
	var buf = new Bytes();
	var f = args.length? new File(args[0]) : new File();
	var tbl = [];
	while (f.readline(buf) >= 0) {
		var t = buf.toString().split("\t");
		var x = parseInt(t[1]);
		if (tbl[x] == null) tbl[x] = 0;
		++tbl[x];
	}
	f.close();
	buf.destroy();
	for (var i = 0; i < tbl.length; ++i)
		if (tbl[i] != null) print(i, tbl[i]);
}

function fm_unreg(args)
{
	var c, min_len = 61, min_utg_len = 0, gap = 1;
	while ((c = getopt(args, 'g:l:u:')) != null) {
		if (c == 'g') gap = parseInt(getopt.arg);
		else if (c == 'l') min_len = parseInt(getopt.arg);
		else if (c == 'u') min_utg_len = parseInt(getopt.arg);
	}

	var f = args.length == getopt.ind? new File() : new File(args[getopt.ind]);
	var b = new Bytes();

	function process(name, len, a)
	{
		var s = 0, e = -1, c = [];
		if (len < min_utg_len) return;
		for (var i = 0; i < a.length; ++i) {
			if (a[i][0] > e) {
				if (e > 0) c.push([s, e]);
				s = a[i][0]; e = a[i][1];
			} else if (a[i][1] > e) e = a[i][1];
		}
		if (e > 0) c.push([s, e]);
		s = 0;
		for (var i = 0; i < c.length; ++i) {
			if (c[i][0] > s && s >= min_len) print(name, s, c[i][0]);
			s = c[i][1];
		}
		if (len > s && len - s >= min_len) print(name, s, len);
		c = [];
	}

	var name, len, a;
	while (f.readline(b) >= 0) {
		var m, line = b.toString();
		if ((m = /^SQ\s(\S+)\s(\d+)/.exec(line)) != null) {
			name = m[1]; len = parseInt(m[2]); a = [];
		} else if ((m = /^EM\s(\d+)\s(\d+)/.exec(line)) != null) {
			var st = parseInt(m[1]), en = parseInt(m[2]);
			if (en - st >= min_len) a.push([st + gap, en - gap]);
		} else if (/^\/\//.test(line)) {
			process(name, len, a);
	}
	}

	b.destroy();
	f.close();
}

/*********************************************
 *** Extract variant sites from raw pileup ***
 *********************************************/

function b8_plp2var(args)
{
	var c, is_fq = false;
	while ((c = getopt(args, 'f')) != null)
		if (c == 'f') is_fq = true;

	var file = args.length > getopt.ind? new File(args[getopt.ind]) : new File();
	var buf = new Bytes();
	var h = { "AC":"M", "CA":"M", "AG":"R", "GA":"R", "AT":"W", "TA":"W", "CG":"S", "GC":"S", "CT":"Y", "TC":"Y", "GT":"K", "TG":"K",
			  "AA":"A", "CC":"C", "GG":"G", "TT":"T" };

	var last_chr = null, last_pos = -1, n_all = 0, sum_all = 0;
	var seq = new Bytes(), qual = new Bytes();
	while (file.readline(buf) >= 0) {
		var t = buf.toString().split("\t");
		if (is_fq && last_chr != t[0]) {
			if (last_chr != null)
				print("@"+last_chr+"\n"+seq.toString()+"\n+\n"+qual.toString()); 
			last_chr = t[0]; last_pos = -1;
		}
		if (!is_fq && !/[ACGTacgt]/.test(t[4])) continue;
		t[2] = t[2].toUpperCase();

		var i = 0, j = 0;
		var alleles = {}, cnt = [], sum = 0;
		while (i < t[4].length && j < t[5].length) {
			var b = t[4].charAt(i);
			if (b == '$') { ++i; continue; }
			if (b == '^') { i += 2; continue; }
			if (b == '*') { ++i, ++j; continue; }

			// determine the allele sequence
			var a, q, forward;
			var match = /[.,A-Za-z]([+-](\d+)[A-Za-z])?/.exec(t[4].substr(i));
			if (b == '.') b = t[2].toUpperCase();
			else if (b == ',') b = t[2].toLowerCase();
			forward = (b.charCodeAt(0) < 97);
			q = t[5].charCodeAt(j) - 33;
			var l_int = 0, l = 0;
			if (match[1] != null) {
				l_int = match[2].length + 1; // including +/-
				l = parseInt(match[2]);
				a = (b + t[4].substr(i + 1, l_int +l)).toUpperCase();
			} else a = b.toUpperCase();
			i += 1 + l_int + l;
			++j;

			// count
			var ci;
			if (alleles[a] == null) alleles[a] = cnt.length, cnt.push([a, 0, 0, 0, 0]);
			ci = alleles[a];
			++cnt[ci][forward? 1 : 2];
			cnt[ci][forward? 3 : 4] += q;
			sum += q;
		}

		sum_all += sum; ++n_all;
		if (!is_fq) {
			var out = [t[0], t[1], t[2], sum, cnt.length];
			for (var i = 0; i < cnt.length; ++i)
				for (var j = 0; j < 5; ++j)
					out.push(cnt[i][j]);
			print(out.join("\t"));
		} else {
			var q = sum <= 93? sum + 33 : 126;
			var pos = parseInt(t[1]) - 1;
			if (last_pos + 1 != pos) {
				for (var i = last_pos + 1; i < pos; ++i) {
					seq.set('N', i); qual.set(33, i);
				}
			}
			if (cnt.length == 1) { // homozygous
				seq.set(cnt[0][0], pos);
			} else if (cnt.length == 2) {
				var a = [cnt[0][0].charAt(0), cnt[1][0].charAt(0)];
				seq.set(h[a[0]+a[1]], pos);
			} else seq.set('N', pos);
			qual.set(q, pos);
			last_pos = pos;
		}
	}
	if (is_fq && last_chr != null) {
		print("@"+last_chr+"\n"+seq.toString()+"\n+\n"+qual.toString());
		warn((sum_all / n_all).toFixed(2));
	}

	buf.destroy();
	file.close();
}

/*************************************
 *** Convert plp2var output to VCF ***
 *************************************/

function b8_var2vcf(args)
{
	var c, qdp = false;
	while ((c = getopt(args, 'q')) != null)
		if (c == 'q') qdp = true;
	var file = args.length > getopt.ind? new File(args[getopt.ind]) : new File();
	var buf = new Bytes();

	while (file.readline(buf) >= 0) {
		var max = 0, match, max_del = '';
		var t = buf.toString().split("\t");
		t[3] = parseInt(t[3]); t[4] = parseInt(t[4]);
		for (var i = 0; i < t[4]; ++i) {
			var match = /^[A-Z]-(\d+)([A-Z]+)/.exec(t[5*(i+1)]);
			if (match != null && max < parseInt(match[1]))
				max = parseInt(match[1]), max_del = match[2];
		}
		var alt = [], dp4 = [0, 0, 0, 0], q = [], qs4 = [0, 0, 0, 0];
		for (var i = 0; i < t[4]; ++i) {
			var a = t[5*(i+1)], match;
			if (a == t[2]) {
				dp4[0] += parseInt(t[5*(i+1) + 1]);
				dp4[1] += parseInt(t[5*(i+1) + 2]);
				qs4[0] += parseInt(t[5*(i+1) + 3]);
				qs4[1] += parseInt(t[5*(i+1) + 4]);
				continue; // identical to the reference
			} else {
				dp4[2] += parseInt(t[5*(i+1) + 1]);
				dp4[3] += parseInt(t[5*(i+1) + 2]);
				qs4[2] += parseInt(t[5*(i+1) + 3]);
				qs4[3] += parseInt(t[5*(i+1) + 4]);
			}
			if ((match = /^[A-Z]\+(\d+)([A-Z]+)/.exec(a)) != null) { // insertion
				alt.push(t[2] + match[2] + max_del);
			} else if ((match = /^[A-Z]-(\d+)([A-Z]+)/.exec(a)) != null) { // deletion
				alt.push(t[2] + max_del.substr(parseInt(match[1])));
			} else { // SNP
				alt.push(a);
			}
			q.push(qs4[2] + qs4[3]);
		}
		if (alt.length == 0) continue; // not a variant
		var alt_sum = 0;
		for (var i = 0; i < q.length; ++i) alt_sum += q[i];
		q.unshift(qs4[0] + qs4[1]);
		var gt;
		if (alt.length == 1 && qs4[0] + qs4[1] == 0) {
			gt = "1/1";
		} else {
			var max = -1, max2 = -1, max_i = -1, max2_i = -1;
			for (var i = 0; i < q.length; ++i) {
				if (max < q[i]) max2 = max, max2_i = max_i, max = q[i], max_i = i;
				else if (max2 < q[i]) max2 = q[i], max2_i = i;
			}
			if (max_i > max2_i) max_i ^= max2_i, max2_i ^= max_i, max_i ^= max2_i;
			gt = max_i + "/" + max2_i;
		}
		var info = qdp? 'DP4='+qs4.join(",")+';RD4='+dp4.join(",") : 'DP4='+dp4.join(",")+';QS4='+qs4.join(",");
		var out = [t[0], t[1], '.', t[2] + max_del, alt.join(","), alt_sum, '.', info, 'GT', gt];
		print(out.join("\t"));
	}

	buf.destroy();
	file.close();
}

/***********************
 *** Main() function ***
 ***********************/

function main(args)
{
	if (args.length == 0) {
		print("\nUsage:    k8 fermi2.js <command> [arguments]\n");
		print("Commands: log2tbl   extract sample info from ropebwt2 log");
		print("          id2sam    relate sequence index to sample info");
		print("          cnt2hist  k-mer count histogram");
		print("          unreg     identify unaligned regions from match output");
		print("          plp2var   extract variant sites from raw pileup");
		print("          var2vcf   convert plp2var output to VCF");
		print("");
		exit(1);
	}

	var cmd = args.shift();
	if (cmd == 'log2tbl') fm_log2tbl(args);
	else if (cmd == 'id2sam') fm_id2sam(args);
	else if (cmd == 'cnt2hist') fm_cnt2hist(args);
	else if (cmd == 'unreg') fm_unreg(args);
	else if (cmd == 'plp2var') b8_plp2var(args);
	else if (cmd == 'var2vcf') b8_var2vcf(args);
	else warn("Unrecognized command");
}

main(arguments);
