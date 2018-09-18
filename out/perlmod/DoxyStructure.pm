sub memberlist($) {
    my $prefix = $_[0];
    return
	[ "hash", $prefix . "s",
	  {
	    members =>
	      [ "list", $prefix . "List",
		[ "hash", $prefix,
		  {
		    kind => [ "string", $prefix . "Kind" ],
		    name => [ "string", $prefix . "Name" ],
		    static => [ "string", $prefix . "Static" ],
		    virtualness => [ "string", $prefix . "Virtualness" ],
		    protection => [ "string", $prefix . "Protection" ],
		    type => [ "string", $prefix . "Type" ],
		    parameters =>
		      [ "list", $prefix . "Params",
			[ "hash", $prefix . "Param",
			  {
			    declaration_name => [ "string", $prefix . "ParamName" ],
			    type => [ "string", $prefix . "ParamType" ],
			  },
			],
		      ],
		    detailed =>
		      [ "hash", $prefix . "Detailed",
			{
			  doc => [ "doc", $prefix . "DetailedDoc" ],
			  return => [ "doc", $prefix . "Return" ],
			  see => [ "doc", $prefix . "See" ],
			  params =>
			    [ "list", $prefix . "PDBlocks",
			      [ "hash", $prefix . "PDBlock",
				{
				  parameters =>
				    [ "list", $prefix . "PDParams",
				      [ "hash", $prefix . "PDParam",
					{
					  name => [ "string", $prefix . "PDParamName" ],
					},
				      ],
				    ],
				  doc => [ "doc", $prefix . "PDDoc" ],
				},
			      ],
			    ],
			},
		      ],
		  },
		],
	      ],
	  },
	];
}

$doxystructure =
    [ "hash", "Root",
      {
	files =>
	  [ "list", "Files",
	    [ "hash", "File",
	      {
		name => [ "string", "FileName" ],
		typedefs => memberlist("FileTypedef"),
		variables => memberlist("FileVariable"),
		functions => memberlist("FileFunction"),
		detailed =>
		  [ "hash", "FileDetailed",
		    {
		      doc => [ "doc", "FileDetailedDoc" ],
		    },
		  ],
	      },
	    ],
	  ],
	pages =>
	  [ "list", "Pages",
	    [ "hash", "Page",
	      {
		name => [ "string", "PageName" ],
		detailed =>
		  [ "hash", "PageDetailed",
		    {
		      doc => [ "doc", "PageDetailedDoc" ],
		    },
		  ],
	      },
	    ],
	  ],
	classes =>
	  [ "list", "Classes",
	    [ "hash", "Class",
	      {
		name => [ "string", "ClassName" ],
		public_typedefs => memberlist("ClassPublicTypedef"),
		public_methods => memberlist("ClassPublicMethod"),
		public_members => memberlist("ClassPublicMember"),
		protected_typedefs => memberlist("ClassProtectedTypedef"),
		protected_methods => memberlist("ClassProtectedMethod"),
		protected_members => memberlist("ClassProtectedMember"),
		private_typedefs => memberlist("ClassPrivateTypedef"),
		private_methods => memberlist("ClassPrivateMethod"),
		private_members => memberlist("ClassPrivateMember"),
		detailed =>
		  [ "hash", "ClassDetailed",
		    {
		      doc => [ "doc", "ClassDetailedDoc" ],
		    },
		  ],
	      },
	    ],
	  ],
	groups =>
	  [ "list", "Groups",
	    [ "hash", "Group",
	      {
		name => [ "string", "GroupName" ],
		title => [ "string", "GroupTitle" ],
		files =>
		  [ "list", "Files",
		    [ "hash", "File",
		      {
		        name => [ "string", "Filename" ]
		      }
		    ],
		  ],
		classes  =>
		  [ "list", "Classes",
		    [ "hash", "Class",
		      {
		        name => [ "string", "Classname" ]
		      }
		    ],
		  ],
		namespaces =>
		  [ "list", "Namespaces",
		    [ "hash", "Namespace",
		      {
		        name => [ "string", "NamespaceName" ]
		      }
		    ],
		  ],
		pages =>
		  [ "list", "Pages",
		    [ "hash", "Page",		      {
		        title => [ "string", "PageName" ]
		      }
		    ],
		  ],
		groups =>
		  [ "list", "Groups",
		    [ "hash", "Group",
		      {
		        title => [ "string", "GroupName" ]
		      }
		    ],
		  ],
		functions => memberlist("GroupFunction"),
		detailed =>
		  [ "hash", "GroupDetailed",
		    {
		      doc => [ "doc", "GroupDetailedDoc" ],
		    },
		  ],
	      }
	    ],
	  ],
      },
    ];

1;
