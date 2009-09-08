// -*- JavaScript -*-
//
// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
//
//                                   Jiao Lin
//                      California Institute of Technology
//                       (C) 2008-2009 All Rights Reserved  
//
// {LicenseText}
//
// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
//


// requires:
//    * luban-core.js


(function(luban, $) {
  
  luban.documentmill = function (actioncompiler) {
    if (actioncompiler == null) {
      actioncompiler = new luban.actioncompiler(this);
    }
    this.actioncompiler = actioncompiler;
  };
  
  luban.documentmill.prototype = {

    'compile': function(action) {
      return this.actioncompiler.compile(action);
    },

    'render': function (doc, parent) {
      if (parent != null) this._parent = parent;
      return this.dispatch(doc);
    },

    'dispatch': function (doc) {
      if (typeof(doc) == 'string') {
	this._parent._je.append(doc);
	return doc;
      }
      var type = doc.type;
      code = 'this.on'+type+"(doc)";
      return eval(code);
    },

    '_onContainer': function(container) {
      var parent = this._parent;
      var type = container.type;
      var factory = luban.elementFactory[type];
      var elem = factory(container, this);
      if (parent != null) parent.add(elem);

      var contents = container.contents;
      if (contents != null) {
	
	for (var i in contents) {
	  this._parent = elem;
	  var subdoc = contents[i];
	  if (typeof(subdoc) == 'string') {
	    elem._je.append(subdoc);
	  } else {
	    subelem = this.dispatch(subdoc);
	  }
	  // elem.add(subelem);
	};

      }
      return elem;
    },
    
    '_onElement': function(element) {
      var type = element.type;
      var factory = luban.elementFactory[type];
      var elem = factory(element, this);
      var parent = this._parent;
      if (parent != null) parent.add(elem);
      return elem;
    },
    
    'onpage': function(page) {
      return this._onContainer(page);
    },

    'oncredential': function(credential) {
      return this._onElement(credential);
    },

    'ondocument': function(document) {
      return this._onContainer(document);
    },

    'onsplitter': function(splitter) {
      var type = splitter.type;
      var factory = luban.elementFactory[type];
      var elem = factory(splitter, this);
      var parent = this._parent;
      if (parent != null) parent.add(elem);

      var contents = splitter.contents;
      if (contents == null) return elem;

      for (var i in contents) {
	this._parent = elem;
	var section = contents[i];
	var subelem = elem.createSection(section);
	
	var subcontents = section.contents;
	if (subcontents == null) continue;
	
	for (var j in subcontents) {
	  this._parent = subelem;
	  var subsubdoc = subcontents[j];
	  this.dispatch(subsubdoc);
	  //subelem.add(this.dispatch(subsubdoc));
	}
      };
      return elem;
    },

    'ontabs': function(tabs) {
      var type = tabs.type;
      var factory = luban.elementFactory[type];
      var elem = factory(tabs, this);
      var parent = this._parent;
      if (parent != null) parent.add(elem);
      elem._je.tabs();

      var contents = tabs.contents
      if (contents == null) return elem;

      var selected_tab=-1;
      for (var i in contents) {
	var tab = contents[i];
	if (selected_tab != -1 && tab.selected) {
	  throw "documentmill.ontabs: multiple tabs are selected";
	}
	if (tab.selected) selected_tab = i;

	var subelem = elem.createTab(tab);
	
	var subcontents = tab.contents;
	if (subcontents == null) continue;
	
	for (var j in subcontents) {
	  var subsubdoc = subcontents[j];
	  subelem.add(this.dispatch(subsubdoc));
	}
      };

      if (selected_tab!=-1) elem._je.tabs('select', selected_tab);

      return elem;
    },

    'ontab': function(tab) {
      elem = this._parent;
      tabelem = elem.createTab(tab);

      var contents = tabelem.contents;
      if (contents != null) {
	
	for (var i in contents) {
	  this._parent = tabelem;
	  var subdoc = contents[i];
	  if (typeof(subdoc) == 'string') {
	    tabelem._je.append(subdoc);
	  } else {
	    this.dispatch(subdoc);
	  }
	};

      }
      return tabelem;
    },

    'onaccordion': function(accordion) {
      var type = accordion.type;
      var factory = luban.elementFactory[type];
      var elem = factory(accordion, this, {plain_element: true});
      var parent = this._parent;
      if (parent != null) parent.add(elem);

      var contents = accordion.contents;
      if (contents == null) return elem;

      var selected_section = null;
      for (var i in contents) {
	var section = contents[i];
	
	this._parent = elem;
	var subelem = this.dispatch(section);
	
	var selected = section.selected;
	if (selected_section!=null && selected) 
	  throw "documentmill.onaccordion: multiple selected sections";
	
	//if (selected) selected_section = section.id;
	if (selected) selected_section = parseInt(i);
	
      };
      
      var opts = $.extend( {}, factory.defaultopts, {active: selected_section});
      elem._je.accordion(opts);
      
      // bind event handler
      var onchange = elem._je.data('onchange-func');
      if (onchange != null)
	elem._je.bind('accordionchange', onchange);
      
      return elem;
    },

    'onaccordionsection': function(section) {
      var parent = this._parent;
      var sectionelem = parent.createSection(section, {plain_element: true});

      var contents = section.contents;
      if (contents != null) {
	
	for (var i in contents) {
	  this._parent = sectionelem;
	  var subdoc = contents[i];
	  if (typeof(subdoc) == 'string') {
	    sectionelem._je.append(subdoc);
	  } else {
	    this.dispatch(subdoc);
	  }
	};

      }
      return sectionelem;
    },

    'onparagraph': function(paragraph) {
      return this._onContainer(paragraph);
    },

    'onlink': function(widget) {
      return this._onElement(widget);
    },

    'ontoolbar': function(toolbar) {
      return this._onContainer(toolbar);
    },

    'ontoolbarspacer': function(widget) {
      return this._onElement(widget);
    },

    'onbutton': function(button) {
      return this._onElement(button);
    },

    'onportlet': function(portlet) {
      return this._onContainer(portlet);
    },

    'onportletitem': function(portletitem) {
      return this._onElement(portletitem);
    },

    'onform': function(form) {
      return this._onContainer(form);
    },

    'onformtextfield': function(formtextfield) {
      return this._onElement(formtextfield);
    },

    'onformpasswordfield': function(formpasswordfield) {
      return this._onElement(formpasswordfield);
    },

    'onformtextarea': function(formtextarea) {
      return this._onElement(formtextarea);
    },

    'onformselectorfield': function(formselectorfield) {
      return this._onElement(formselectorfield);
    },

    'onformsubmitbutton': function(formsubmitbutton) {
      return this._onElement(formsubmitbutton);
    },

    'onappmenubar': function(appmenubar) {
      var type = appmenubar.type;
      var factory = luban.elementFactory[type];
      var elem = factory(appmenubar, this);
      var parent = this._parent;
      if (parent != null) parent.add(elem);

      var contents = appmenubar.contents
      if (contents == null) return elem;

      for (var i in contents) {
	this._parent = elem;
	var subdoc = contents[i];
	this.dispatch(subdoc);
      };

      elem._je.jdMenu();
      return elem;
    },

    'onappmenu': function(widget) {
      return this._onContainer(widget);
    },

    'onappmenuitem': function(widget) {
      return this._onElement(widget);
    },

    'ontreeview': function(treeview) {
      var type = treeview.type;
      var factory = luban.elementFactory[type];
      var elem = factory(treeview, this);
      var parent = this._parent;
      if (parent != null) parent.add(elem);

      var root = treeview.root;
      if (root) {
	this._parent = elem;
	this.dispatch(root);
      }

      var rules = {};
      if (treeview.draggable)
	rules.draggable = 'all';

      var callback = {};
      if (treeview.onnodemoving) 
	callback.onmove = elem._je.data('onnodemoving-func');

      opts = {
	  'rules': rules,
	  'callback': callback,
	}
      elem._je.tree(opts);
      // $.tree_reference(treeview.id).open_all();
      // elem._je.treeview();
      return elem;
    },

    'ontreeviewbranch': function(widget) {
      return this._onContainer(widget);
    },

    'ontreeviewleaf': function(widget) {
      return this._onElement(widget);
    },

    'ontable': function(widget) {
      var type = widget.type;
      var factory = luban.elementFactory[type];
      var elem = factory(widget, this._parent, this);
      return elem;
    }

  };


  luban.docmill = new luban.documentmill;

 })(luban, jQuery);


// End of file
