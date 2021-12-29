import sqlite3
from numbers import Number
from functools import reduce

class SplitFilter():
    TAX_KEYWORDS = ['SUPERKINGDOM', 'PHYLUM', 'CLASS', 'ORDER', 'FAMILY', 'GENUS', 'SPECIES']
    FUN_KEYWORDS = ['FUN', 'FUNH']
    SAMPLES = []
    ALL_KEYWORDS = TAX_KEYWORDS + FUN_KEYWORDS + SAMPLES
    RELATIONAL_OPERATORS = ['==', '!=', '>=', '<=', '>', '<', 'IN', 'NOT IN', 'CONTAINS', 'DOES NOT CONTAIN']
    LOGICAL_OPERATORS = ['AND', 'OR']
    ALL_OPERATORS = RELATIONAL_OPERATORS + LOGICAL_OPERATORS

    # NEW
    def __init__(self, contigs_db_path, profile_db_path, contigs_tax_path, anvio_version):
    #def __init__(self, contigs_db_path, profile_db_path, contigs_tax_path):
    #
        """
        Parse a string cointaining relational expressions joined by logical operators.
        """
        self.contigs_db_path  = contigs_db_path
        self.profile_db_path  = profile_db_path
        self.contigs_tax_path = contigs_tax_path
        # NEW
        self.anvio_version = anvio_version
        #
        self.has_blank_profile = self.is_blank_profile(self.profile_db_path)
        if not self.has_blank_profile:
            self.SAMPLES = self.load_anvio_samples(self.profile_db_path)
            # NEW
            self.splits  = self.load_anvio_splits(self.profile_db_path, 'abundance_splits', anvio_version)
            #self.splits  = self.load_anvio_splits(self.profile_db_path, 'abundance_splits')
            #
        else:
            self.SAMPLES = []
            # NEW
            self.splits  = self.load_anvio_splits(self.profile_db_path, 'atomic_data_splits', anvio_version)
            #self.splits  = self.load_anvio_splits(self.profile_db_path, 'atomic_splits')
            #
        self.ALL_KEYWORDS = self.TAX_KEYWORDS + self.FUN_KEYWORDS + self.SAMPLES # Add sample names to the list of valid keywords.
        self.contigTax = self.load_taxonomy(self.contigs_tax_path)


    @staticmethod
    def is_blank_profile(profile_db_path):
        """
        From our tests, the profile database contains the table "abundance_splits"
        the profile has samples, and the table "atomic_data_splits" if the profile
        is blank.
        """ 
        conn = sqlite3.connect(profile_db_path)
        c = conn.cursor()
        c.execute('SELECT * FROM SQLITE_MASTER')
        tables = [x[2] for x in c.fetchall() if x[0]=='table']
        if 'atomic_data_splits' in tables and 'abundance_splits' in tables:
            raise Exception('Your profile database contains one table named "atomic_data_splits" and another named "abundance_splits", \
                             so we can\'t tell whether it is a blank profile or not. This might be our fault, let us know.')
        return 'atomic_data_splits' in tables
    


    @staticmethod
    def load_anvio_samples(profile_db_path):
        conn = sqlite3.connect(profile_db_path)
        c = conn.cursor()
        #c.execute('PRAGMA table_info("abundance_splits")')
        #samples = [x[1] for x in c.fetchall()][1:-1]
        # NEW
        # Affected by the new structure of PROFILE.db in anvio 7.1. Read from self?
        c.execute('SELECT * FROM self WHERE key == \'samples\'')
        samples = list(map(list, c.fetchall()))[0]
        conn.close()
        return samples

    

    @staticmethod
    # NEW
    def load_anvio_splits(profile_db_path, table, anvio_version):
    #def load_anvio_splits(profile_db_path, table):
    #
        conn = sqlite3.connect(profile_db_path)
        c = conn.cursor()
        c.execute('SELECT * FROM SQLITE_MASTER')
        # NEW
        if anvio_version > 7:
            c.execute('SELECT item FROM {}'.format(table)) # named contig, but they're really splits
        else:
            c.execute('SELECT contig FROM {}'.format(table)) # named contig, but they're really splits
        #c.execute('SELECT contig FROM {}'.format(table))
        #
        splits = set([x[0] for x in c.fetchall()])
        conn.close()   
        return splits



    @staticmethod
    def load_taxonomy(contigs_tax_path):
        contigTax = {}
        with open(contigs_tax_path) as infile:
            infile.readline()
            for line in infile:
                line = line.strip().split('\t')
                assert len(line) == 8
                contigTax[line[0]] = line[1:]
        return contigTax



    def get_parsed_annotations(self):
        return self.splits, self.SAMPLES, self.contigTax
    


    def parse_query(self, query):
        """
        Generate a syntax tree from a query expression.
        
        From:
            (PHYLUM IN [Bacteroidetes] AND CLASS!=Flavobacteria) OR CLASS IN [Alphaproteobacteria, Gammaproteobacteria] OR (HADZA3_SRR1927149>0.5 AND HADZA3_SRR1927149 < 1)
        To:
             OR
                 AND
                         ('IN', 'PHYLUM', '[Bacteroidetes]')
                         ('!=', 'CLASS', 'Flavobacteria')
                 ('IN', 'CLASS', '[Alphaproteobacteria, Gammaproteobacteria]')
                 AND
                         ('>', 'HADZA3_SRR1927149', '0.5')
                         ('<', 'HADZA3_SRR1927149', '1')
        """
        return self.parse_logical(self.parse_parentheses(query))


    @staticmethod
    def parse_parentheses(input_string):
        """
        Parse a raw string into a list of subexpressions that were originally separated by parentheses.

        From:
            (PHYLUM IN [Bacteroidetes] AND CLASS!=Flavobacteria) OR CLASS IN [Alphaproteobacteria, Gammaproteobacteria] OR (HADZA3_SRR1927149>0.5 AND HADZA3_SRR1927149 < 1)
        To:
            ['PHYLUM IN [Bacteroidetes] $AND$ CLASS!=Flavobacteria', '$OR$ CLASS IN [Alphaproteobacteria, Gammaproteobacteria] $OR$', 'HADZA3_SRR1927149>0.5 $AND$ HADZA3_SRR1927149 < 1']
        """
        tree = {'root': {'parent': None, 'children': []}}
        current = 'root'
        buff=''
        i = 0

        ST = input_string.replace(' OR ', ' $OR$ ').replace(' AND ', ' $AND$ ').replace(' OR(', ' $OR$ (').replace(' AND(', ' $AND$ (')
        ST = ST.replace('IN[', 'IN [')
        pars = 0
        for x in ST:
            if x == '(':
                pars += 1
                i+=1
                buff = buff.strip()
                if buff:
                    tree[current]['children'].append(buff)
                    buff = ''
                new = []
                tree[current]['children'].append(new)
                tree[i] = {'parent': current, 'children': new}
                current = i
            elif x==')':
                pars -= 1
                buff = buff.strip()
                if buff:
                    tree[current]['children'].append(buff)
                    buff = ''
                # If a clause has only one item, just store the string instead of a one-item list.
                if len(tree[current]['children'])==1:
                    tree[tree[current]['parent']]['children'][-1] = tree[current]['children'][0]
                current = tree[current]['parent']
            else:
                buff += x
        buff = buff.strip()
        if buff:
            tree[current]['children'].append(buff)

        if pars > 0:
            raise Exception('Your expression is missing a closing parenthesis (or has too many opening parentheses!)')
        if pars < 0:
            raise Exception('Your expression is missing an opening parenthesis (or has too many closing parentheses!)')

        return tree['root']['children']
    

    def parse_logical(self, clause):
        """
        Parse the subexpressions generated by parse_parentheses into a syntax tree.
        
        """
        
        def parse_relational(expr):
            """
            Parse a string consisting a single relational expression into a tuple of the form (operation, subject, value)
            """
            if expr:
                temp_expr = expr.replace(' NOT IN ', ' $not in$ ').replace(' CONTAINS ', ' $contains$ ').replace(' DOES NOT CONTAIN ', ' $does not contain$ ').replace(' IN ', ' $in$ ')
                # Hack since "IN" is itself a substring of "CONTAINS", "SUPERKINGDOM",
                #  "NOT IN" and "DOES NOT CONTAIN"
                temp_relOps = [op.replace('NOT IN', '$not in$').replace('CONTAINS', '$contains$').replace('DOES NOT CONTAIN', '$does not contain$').replace('IN', '$in$')
                               for op in self.RELATIONAL_OPERATORS] # Keep hacking
                ops = [op.replace('$','').upper() for op in temp_relOps if op in temp_expr] # Finish hacking
                if len(ops) != 1:
                    raise Exception('Either none or more than one relational operators in expression "{}"'.format(expr))
                op = ops[0]
                subject, value = [x.strip() for x in expr.split(op)]
                if subject not in self.ALL_KEYWORDS:
                    raise Exception('Unknown keyword "{}" in expression "{}"'.format(subject, expr))
                assert subject in self.ALL_KEYWORDS
                if op in ('IN', 'NOT IN'):
                    if not value.startswith('[') or not value.endswith(']'):
                        raise Exception('Syntax error in expression "{}"'.format(expr))
                    value = set([x.strip() for x in value[1:-1].split(',')])
                return (op, subject, value)
            else:
                return None

        if type(clause) is str:
            isOR = '$OR$' in clause
            isAND = '$AND$' in clause
            isLogic = isOR or isAND
            if isOR and isAND:
                raise Exception('Clause: "{}" can not contain both AND and OR'.format(clause))
            if isOR:
                op, split = 'OR', '$OR$'
            elif isAND:
                op, split = 'AND', '$AND$'
            if isLogic:
                parsed = [op] + [[parse_relational(x) for x in clause.split(split)]]
            else:
                parsed = parse_relational(clause)

        elif type(clause) is list:
            # The parse parentheses function will separate expressions joined by the same logical operator into two tuple elements.
            # eg from this string
            #   '(PHYLUM IN [Bacteroidetes] AND CLASS!=Flavobacteria) OR CLASS IN [Alphaproteobacteria, Gammaproteobacteria] OR (HADZA3_SRR1927149>0.5 AND HADZA3_SRR1927149 < 1)'
            # it generates:
            #   ['PHYLUM IN [Bacteroidetes] $AND$ CLASS!=Flavobacteria', '$OR$ CLASS IN [Alphaproteobacteria, Gammaproteobacteria] $OR$', 'HADZA3_SRR1927149>0.5 $AND$ HADZA3_SRR1927149 < 1']
            # note how the middle subexpression starts and ends with $OR$.
            # second_parts and first_parts detect whether a subexpression is actually the second or first part of a longer logical expression, and generates the right tree.
            parsed = []
            second_parts = [i for i, sub in enumerate(clause) if type(sub) is str and (sub.startswith('$OR$') or sub.startswith('$AND$'))]
            first_parts = [i for i, sub in enumerate(clause) if type(sub) is str and (sub.endswith('$OR$') or sub.endswith('$AND$'))]
            for i, sub in enumerate(clause):
                if i-1 in first_parts or i+1 in second_parts:
                    continue
                if i in first_parts and i in second_parts:
                    op = 'OR' if sub.startswith('$OR$') else 'AND'
                    op2 = 'OR' if sub.endswith('$OR$') else 'AND'
                    split = '${}$'.format(op)
                    if op!=op2:
                        raise Exception('Clause: "{}" can not contain both AND and OR'.format(sub))
                    if sub.strip() == split:
                        parsed.append([op, [self.parse_logical(clause[i-1]), self.parse_logical(clause[i+1])]])
                    else:
                        sub = sub.split(split,1)[1].rsplit(split,1)[0]
                        parsed.append([op, [self.parse_logical(clause[i-1]), self.parse_logical(sub), self.parse_logical(clause[i+1])]])
                elif i in first_parts:
                    op = 'OR' if sub.endswith('$OR$') else 'AND'
                    split = '${}$'.format(op)
                    sub = sub.rsplit(split,1)[0]
                    parsed.append([op, [self.parse_logical(sub), self.parse_logical(clause[i+1])]])
                elif i in second_parts:
                    op = 'OR' if sub.startswith('$OR$') else 'AND'
                    split = '${}$'.format(op)
                    sub = sub.split(split,1)[1]
                    parsed.append([op, [self.parse_logical(clause[i-1]), self.parse_logical(sub)]])
                else:
                    parsed.append(self.parse_logical(sub))
                
        if len(parsed) == 1:
            parsed = parsed[0]
        return parsed
        

    def print_tree(self, query):
        tree = self.parse_query(query)
        self._print_node(tree)


    def _print_node(self, node, level=0):
        """
        A node can have the following structures:
        - [LOGICAL_OP, [list of expressions]] -> evaluate each expression, apply set of splits that fulfill LOGICAL_OP on all of them.
        - (RELATIONAL_OP, subject, value)     -> evaluate relational op and return set of splits which fulfill it.
        """
        if node[0] in self.LOGICAL_OPERATORS:
            print('{}{}'.format('\t'*level, node[0]))
            for sub in node[1]:
                self._print_node(sub, level+1)
        elif node[0] in self.RELATIONAL_OPERATORS:
            print('{}{}'.format('\t'*level, node))
        else:
            raise Exception('This should not happen')


    def run(self, query):
        tree = self.parse_query(query)
        goodSplits = self._run_node(tree)
        goodSplits = goodSplits & self.splits # Select only splits that are in the profile database
        return goodSplits


    def _run_node(self, node):
        """
        A node can have the following structures:
        - [LOGICAL_OP, [list of expressions]] -> evaluate each expression, apply set of splits that fulfill LOGICAL_OP on all of them.
        - (RELATIONAL_OP, subject, value)     -> evaluate relational op and return set of splits which fulfill it.
        """
        if node[0] in self.LOGICAL_OPERATORS:
            goodSplits = map(self._run_node, node[1])
            if node[0] == 'AND':         
                return reduce(lambda x,y: x & y, goodSplits)
            elif node[0] == 'OR':
                return reduce(lambda x,y: x | y, goodSplits)
        elif node[0] in self.RELATIONAL_OPERATORS:
            return self.runExpr(*node)
        else:
            raise Exception('This should not happen')


    def runExpr(self, op, subject, value):
        if subject in self.TAX_KEYWORDS:
            res = self.runTax(op, subject, value)
        elif subject in self.FUN_KEYWORDS:
            res = self.runFun(op, subject, value)
        elif subject in self.SAMPLES:
            res = self.runAbund(op, subject, value)
        return res


    def runTax(self, op, subject, value):
        assert op not in ('>', '<', '>=', '<=')
        taxIdx = {rank: i for i, rank in enumerate(self.TAX_KEYWORDS)}
        goodTax = set()
        for split in self.splits:
            contig = split.rsplit('_split', 1)[0]
            tax = self.contigTax[contig][taxIdx[subject]]
            if self.evaluate(op, tax, value):
                goodTax.add(split)
        return goodTax


    def runFun(self, op, subject, value):
        assert op in ('==', '!=', 'IN', 'NOT IN', 'CONTAINS', 'DOES NOT CONTAIN')
        conn = sqlite3.connect(self.contigs_db_path)
        c = conn.cursor()
        if op in ('IN', 'NOT IN'):
            assert type(value) is set
            valueSub = '({})'.format(','.join('?' * len(value) ))
            value = 2 * list(value)
        else:
            valueSub = '?'
            if op == 'CONTAINS':
                op = 'LIKE'
                value = '%{}%'.format(value)
            if op == 'DOES NOT CONTAIN':
                op = 'NOT LIKE'
                value = '%{}%'.format(value)
            value = [value, value]
        if subject == 'FUN':
            sources = '"KEGG", "COG", "PFAM"'
        elif subject == 'FUNH':
            sources = '"KEGGPATH"'
        else:
            raise Exception('This should not happen')
        query = 'SELECT DISTINCT genes_in_splits.split FROM genes_in_splits, gene_functions WHERE gene_functions.source IN ({}) AND (gene_functions.function {} {} OR gene_functions.accession {} {}) AND gene_functions.gene_callers_id==genes_in_splits.gene_callers_id;'.format(sources, op, valueSub, op, valueSub)
        print(query)
        # DISTINCT is not really needed as we will cast results into a set, but is clearer this way.
        c.execute(query, value)
        goodSplits = set([x[0] for x in c.fetchall()])
        conn.close()
        return goodSplits
        

    def runAbund(self, op, subject, value):
        assert op in ('==', '!=', '>', '>=', '<', '<=')
        assert subject in self.SAMPLES
        try:
            float(value)
        except:
            raise Exception('Value in "{} {} {}" must be numeric!'.format(subject, op, value))
        conn = sqlite3.connect(self.profile_db_path)
        c = conn.cursor()
        #op, subject and value should have just been tested for sanity so we do direct substitution.
        c.execute('SELECT contig FROM abundance_splits WHERE {}{}{};'.format(subject, op, value)) # we say SELECT contig but they're really splits
        goodSplits = set([x[0] for x in c.fetchall()])
        conn.close()
        return goodSplits


    @staticmethod
    def evaluate(op, tSubject, tValue):
        #tSubject and tValue are true subjects and values, directly comparable after having been parsed and casted to the appropriate types.
        if op == '==':
            return tSubject == tValue
        elif op == '!=':
            return tSubject != tValue
        elif op == 'IN':
            return tSubject in tValue
        elif op == 'NOT IN':
            return tSubject not in tValue
        elif op == 'CONTAINS':
            return tValue in tSubject
        elif op == 'DOES NOT CONTAIN':
            return tValue not in tSubject
        elif op == '>':
            return tSubject > tValue
        elif op == '>=':
            return tSubject >= tValue
        elif op == '<':
            return tSubject < tValue
        elif op == '<=':
            return tSubject <= tValue
        else:
            raise Exception('Unrecognized operation "{}"'.format(op))



class DictFilter(SplitFilter):
    def __init__(self, read_tax, tax_source, read_fun, fun_info, tax2fun, fun2tax):
        self.read_tax = read_tax[tax_source]
        self.read_fun = read_fun
        self.fun_info = fun_info
        self.tax2fun  = tax2fun
        self.fun2tax  = fun2tax
        self.ALL_KEYWORDS = self.TAX_KEYWORDS + self.FUN_KEYWORDS


    def run(self, query):
        tree = self.parse_query(query)
        goodReads = self._run_node(tree)
        return goodReads


    def runTax(self, op, subject, value):
        assert op not in ('>', '<', '>=', '<=')
        taxIdx = {rank: i for i, rank in enumerate(self.TAX_KEYWORDS)}
        goodTax = set()
        for read, tax in self.read_tax.items():
            tax = tax[taxIdx[subject]].split(';')[-1].split('_', 1)[1] # get rid of "k_" in "k_Bacteria"
            if self.evaluate(op, tax, value):
                goodTax.add(read)
                goodTax.update(self.tax2fun[read])
        return goodTax


    def runFun(self, op, subject, value):
        assert op in ('==', '!=', 'IN', 'NOT IN', 'CONTAINS', 'DOES NOT CONTAIN')
        goodFun = set()
        hmethods = {'kegg', 'cogs'}
        for method in self.read_fun:
            if subject == 'FUNH' and method not in hmethods:
                continue
            for read, funs in self.read_fun[method].items():
                for fun_id in funs:
                    # Get info
                    if method in hmethods:
                        if fun_id in self.fun_info[method]: 
                            fun_name, fun_h = self.fun_info[method][fun_id]
                        else:
                            fun_name, fun_h = None, None
                    else:
                        fun_name = self.fun_info[method][fun_id] if fun_id in self.fun_info[method] else None
                    # Select relevant info
                    if subject == 'FUN':
                        info = (fun_id, fun_name) if fun_name else (fun_id,)
                    elif subject == 'FUNH':
                        info = (fun_h,) if fun_h else tuple()
                    # Filter
                    for i in info:
                        if self.evaluate(op, i, value):
                            goodFun.add(read)
                            goodFun.add(self.fun2tax[read])
        return goodFun
                   

            
