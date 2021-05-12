## Software installation for storage and reporting frontend

### Requirements
* [Apache 2.4](https://httpd.apache.org)
  * must be cgi enabled (serve-cgi-bin.conf)

Some other webserver could be used but the file organization and
security model for CLAMP were design to fit into a pre-existing
server running a variety of applications. An example httpd.conf
follows (where "path/to..." should obviously be localized).

```
<Directory path/to/cml>
AuthName Foo
AuthType Basic
AuthUserFile path/to/passwd
AuthGroupFile path/to/passwd/group
AddHandler cgi-script .cgi
Options +Indexes +ExecCGI +FollowSymLinks
DirectoryIndex index.cgi
AllowOverride Limit
Require valid-user
</Directory>

<Directory path/to/cml/reg>
Require group cml
</Directory>
```

The _password_ file looks like:
```
alice:<encrypted>
bob:<encrypted>
carol:<encrypted>
```

The _group_ file looks like:
```
cml: bob carol
```

So, in this example, Bob and Carol would be able to register and upload samples.
Alice would only be allowed to search and view the results.


* [Perl](https://www.perl.org), eg version 5.22 but at least 5.10
  * Perl Modules (with minimum or suggested version):
  * CGI 4.26
  * CGI::Carp 4.26
  * Compress::Zlib 2.068
  * DBD::SQLite 1.56
  * File::Copy 2.30
  * GD 2.53
  * GD::Graph 1.54
  * GD::Text 0.86
  * List::Util 1.41
  * Time::Piece 1.29
  * Time::Seconds 1.29

Note that some of these Perl modules are or have been in Perl core.
Some modules may have dependencies (most notably _GD_).
The given version numbers are known to work but are unlikely to be 
strictly required.
