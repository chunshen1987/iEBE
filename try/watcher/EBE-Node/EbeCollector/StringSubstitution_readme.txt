
==============================================
    Document for StringSubstitution Module
==============================================

The StringSubstitution module provides a StringSubstitution class to support the idea of performming string substituion using rules, which is explained in the following. The StringSubstitution class uses one dictionary to specify all the replacement rules, which can be set either using constructor or by an explicit call to the setRules member function. A single replacement rule can be repeatedly applied by using the applySingleRule function, which, together with other functions, will be explained in details in the following.


1) "Rule" for string substitution.

A "rule" is a ("key", "value") pair. A "key" can be either a string or a compiled re pattern (faster in performance), and a "value" can be either a string or a function, and it is used to generate the string used for substituion.

If the "value" is a string, then it is formatted against with groups() list of the resulting matched object to generate the substitution string. If the "value" is a function, then it is called with a single argument, the matched object, to generate the substitution string. If it is neither a string or a function, it will be converted to a string using the str() function.

The generated substitution string will replace group(0) of the matched object for its occurance in the target string; such an operation is called a string substitution using a rule.


2) Apply a single rule.

The static applySingleRule(rule, expression, restrictNumberOfRecursionsTo=None) function applies a single "rule". This function apply the substitution recursively and return the resulting string and the number of recursions of substitutions.

The following are a few examples.
>>> from StringSubstitution import StringSubstitution

Simple replacement of all letter "s" by "w" in a string:
>>> StringSubstitution.applySingleRule(("s","w"), "A string is a string.")
('A wtring iw a wtring.', 1)

Using a formatted string as "value": search for "s" with the closest "r", then replace them by "w" and "u" pair using groups. Note that in order for the formatting to work, the target string needs to contain markers like "{0[#]}":
>>> StringSubstitution.applySingleRule(("s(.*?)r", "w{0[0]}u"), "sxr sxxr sxxxr")
('wxu wxxu wxxxu', 3)

Using a function as "value": adding the two adjacent single digits into one (1234->334->64->10->1):
>>> StringSubstitution.applySingleRule(("(\d).*?(\d)", lambda m: "%d" % (int(m.group(1))+int(m.group(2)))), "1234")
('1', 4)

When restrictNumberOfRecursionsTo is set to an integer, only the specified number of recursion of substitutions will be performed. The following example repeat the example above, but perform only 1 recursion of substituion:
>>> StringSubstitution.applySingleRule(("(\d).*?(\d)", lambda m: "%d" % (int(m.group(1))+int(m.group(2)))), "1234", restrictNumberOfRecursionsTo=1)
('334', 1)


3) Apply a set of rules.

The main purpose of the StringSubstitution class is to store a set of compiled rules so that it can be applied efficiently to any given expressions. The set of rules usually is a lengthy list of (key, value) tuple and it should be first provided to the class, which will then be compiled for efficiency. To supply a set of rules either use the constructor or use the setRules function.

The following is an example. Suppose we want to change all "a" to "1", "b" to "2", "c" to "3", then we can initialize a class using the following:
>>> substitutor = StringSubstitution( (("a",1), ("b",2), ("c",3)) )

Then the applyAllRules(expression) function can be used formulate a string using all the rules. For example to formulate the string "abcaabbccccccc" (for those who know what this means: good old time...):
>>> substitutor.applyAllRules("abcaabbccccccc")
('12311223333333', 1)

The rules can be change by the setRules(rules) function, for example:
>>> substitutor.setRules( (("1","a"), ("2","b"), ("3","c")) )
>>> substitutor.applyAllRules("12311223333333")
('abcaabbccccccc', 1)

Note that the class works in the boundary case when there is no rules given:
>>> StringSubstitution().applyAllRules("Hello StringSubstitution!")
('Hello StringSubstitution!', 0)

The applyAllRules function can accept two more arguments restrictNumberOfRecursionsPerScan and restrictNumberOfScans. The former is used as the restrictNumberOfRecursionsTo argument when calling the applySingleRule function when applying all rules; and the later can be used to restrict the total number of loop scans all rules are applied. As a dumb example, suppose we have two rules, one changing "a" to "b" and one changing "b" to "a", then without restricting the number of loop scans the substitution will never end:
>>> StringSubstitution( (("a","b"), ("b","a")) ).applyAllRules("a", restrictNumberOfScans=99)
('a', 99)



===========
The END
===========
