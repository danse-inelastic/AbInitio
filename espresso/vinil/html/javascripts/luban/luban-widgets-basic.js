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


  // aliases
  var ef = luban.elementFactory;
  var widgets = luban.widgets;


  // page
  ef.page = function (kwds) {
    var title = kwds.title;
    document.title = title;

    var ret = $('#body\-wrapper');
    var id = kwds.id;
    if (id) ret.attr('id', id);
    ret.addClass('luban-page');

    ret.data('luban-element-type', 'page');
    
    return ret.lubanElement();
  };
  
  widgets.page = function(elem) {
    this.super1 = widgets.base;
    this.super1(elem);
  };
  widgets.page.prototype = new widgets.base;
  widgets.page.prototype.setAttribute = function(attrs) {
      var title = attrs.title;
      if (title) this.setTitle(title);
  };
  widgets.page.prototype.setTitle = function (title) {
    document.title = title;
  };
  widgets.page.prototype.empty = function() {
    this._je.empty();
  };


  // document
  ef.document = function (kwds, docmill) {
    var Class = kwds.Class;
    var id = kwds.id;
    var ret = tag('div', {'id': id});
    ret.addClass(Class);

    var titlediv = tag('h1');
    titlediv.addClass(widgets.document.prototype.title_class);
    ret.append(titlediv);

    var title = kwds.title;
    if (title) 
      titlediv.text(title);
    else
      titlediv.hide();

    var onclick = kwds.onclick;
    if (onclick) 
      ret.click( function() { docmill.compile(onclick); return false; } );

    return ret.lubanElement('document');
  };

  widgets.document = function(elem) {
    this.super1 = widgets.base;
    this.super1(elem);
  };
  widgets.document.prototype = new widgets.base;
  widgets.document.prototype.title_class = 'luban-document-title';
  widgets.document.prototype.empty = function() {
    // keep the title
    var div = this._je;
    var titlediv = div.children('.'+this.title_class);
    this._je.empty();
    if (titlediv != null)
      this._je.append(titlediv);
  };
  widgets.document.prototype.setAttribute = function(attrs) {
    var div = this._je;

    var title = attrs.title;
    var titlediv = div.children('.'+this.title_class);
    if (title) {
      titlediv.text(title).show();
    } else {
      if (title=='') 
	titlediv.text(title).hide();
    }

    var Class = attrs.Class;
    if (Class) {
      div.removeClass();
      div.addClass(Class);
    }
  };

  // splitter
  //  factory
  ef.splitter = function(kwds, docmill) {
    var orientation = kwds.orientation;
    var id = kwds.id;
    
    var container = tag('div', {"id": id});
    var kls = kwds.Class;
    if (kls) container.addClass(kls);

    if (orientation == 'vertical') {
    } else {
      table = tag('table'); container.append(table);
      table.addClass('luban-splitter');
      tbody = tag('tbody'); table.append(tbody);
      tr = tag('tr'); tr.addClass('luban-splitter-section-container');
      tbody.append(tr);
    }

    var onclick = kwds.onclick;
    if (onclick) 
      container.click( function() { docmill.compile(onclick); return false; } );

    container.data('orientation', orientation);
    return container.lubanElement('splitter');
  };

  //  object
  widgets.splitter = function(elem) {
    this.super1 = widgets.base;
    this.super1(elem);
  };
  widgets.splitter.prototype = new widgets.base;  
  widgets.splitter.prototype.createSection = function(kwds, docmill) {
    var id = kwds.id;
    var kls = kwds.Class;
    var onclick = kwds.onclick;

    var section_container_class = 'luban-splitter-section-container';

    var orientation = this._je.data('orientation');
    var tagname, sec, interior_container;
    if (orientation == 'vertical') {
      tagname = 'div';
      sec = tag(tagname, {'id': id});
      interior_container = this._je;
    } else {
      tagname = 'td';
      sec = tag(tagname, {'id': id});
      // sec.addClass('luban-splitter-section'); why need this? already got a classifier in the end at sec.lubanElement('splitsection')
      var table = this._je.children('table');
      var tbody = table.children('tbody');
      interior_container = tbody.children('tr.'+section_container_class);
    }
    
    interior_container.append(sec);

    var div = tag('div');
    div.addClass(widgets.splitsection.prototype.interior_container_class);
    sec.append(div);
    
    if (kls) div.addClass(kls);

    if (onclick)
      div.click( function() { docmill.compile(onclick); return false; } );

    return sec.lubanElement('splitsection');
  };

  // splitsection
  widgets.splitsection = function(elem) {
    this.super1 = widgets.base;
    this.super1(elem);
  };
  widgets.splitsection.prototype = new widgets.base;
  widgets.splitsection.prototype.interior_container_class = 'luban-splitsection-interior-container';
  widgets.splitsection.prototype.add = function(elem) {
    var c = this.jqueryelem.children('div.'+this.interior_container_class);
    c.append(elem.jqueryelem);
  };
  widgets.splitsection.prototype.addClass = function(Class) {
    var c = this.jqueryelem.children('div'+this.interior_container_class);
    c.addClass(Class);
  };
  widgets.splitsection.prototype.removeClass = function(Class) {
    var c = this.jqueryelem.children('div'+this.interior_container_class);
    if (c.hasClass(Class)) {
	c.removeClass(Class);
    }
  };

  // paragraph
  //  factory
  ef.paragraph = function(kwds, docmill) {
    var id = kwds.id;
    var p = tag('p', {"id":id});
    p.text(kwds.text.join('\n'));

    var kls = kwds.Class;
    if (kls) p.addClass(kls);

    var onclick = kwds.onclick;
    if (onclick) 
      p.click( function() { docmill.compile(onclick); return false; } );

    return p.lubanElement('paragraph');
  };
  //  object
  widgets.paragraph = function(elem) {
    this.super1 = widgets.base;
    this.super1(elem);
  };
  widgets.paragraph.prototype = new widgets.base;


  // htmldocument
  //  factory
  ef.htmldocument = function(kwds, docmill) {
    var id = kwds.id;
    var div = tag('div', {"id":id});
    div.html(kwds.text.join('\n'));

    var kls = kwds.Class;
    if (kls) div.addClass(kls);

    var onclick = kwds.onclick;
    if (onclick) 
      div.click( function() { docmill.compile(onclick); return false; } );

    return div.lubanElement('htmldocument');
  };
  //  object
  widgets.htmldocument = function(elem) {
    this.super1 = widgets.base;
    this.super1(elem);
  };
  widgets.htmldocument.prototype = new widgets.base;

  // link
  //  factory
  ef.link = function(kwds, docmill) {
    var id = kwds.id;
    var a = tag('a', {"id":id});
    a.append(kwds.label);
    a.addClass('luban-link');
    
    var kls = kwds.Class;
    if (kls) a.addClass(kls);

    var onclick = kwds.onclick;
    if (onclick != null && onclick != '') 
      a.click( function() { docmill.compile(onclick); return false; } );
    return a.lubanElement('link');
  };
  //  object
  widgets.link = function(elem) {
    this.super1 = widgets.base;
    this.super1(elem);
  };
  widgets.link.prototype = new widgets.base;


  // button
  //  factory
  ef.button = function(kwds, docmill) {
    var id = kwds.id;
    var div = tag('div', {"id": id}); div.addClass('luban-button-container');
//     var table = tag('table');  table.addClass('luban-button-container');
//     div.append(table);
//     var tr = tag('tr'); table.append(tr);
//     var td = tag('td'); td.addClass('luban-button'); tr.append(td);
    var a = tag('a', {"title": kwds.tip}); //td.append(a);
    div.append(a); a.addClass('luban-button');

    var icon = kwds.icon; var img;
    if (icon != null && icon != '') {
      img = tag('img', {'src': luban.iconpath(icon)});
      a.append(img);
      img.addClass('luban-buttonIcon');
    }

    var label = kwds.label; var span;
    if (label != null && label != '') {
      span = tag('span');
      a.append(span);
      span.text(label);
    }

    var onclick = kwds.onclick;
    if (onclick != null && onclick != '')
      a.click( function() { docmill.compile(onclick); return false; } );

    div.data('luban-element-type', 'button');
    return div.lubanElement();
  };
  //  object
  widgets.button = function(elem) {
    this.super1 = widgets.base;
    this.super1(elem);
  };
  widgets.button.prototype = new widgets.base;



  // progressbar
  //  factory
  ef.progressbar = function(kwds, docmill) {
    var id = kwds.id;
    var div = tag('div', {"id": id}); div.addClass('luban-progressbar');

    var kls = kwds.Class;
    if (kls) div.addClass(kls);

    // periodically trigger a "checking" event
    var onchecking = kwds.onchecking;
    if (!onchecking) throw "need to define onchecking handler";

    var onfinished = kwds.onfinished;
    if (!onfinished) throw 'need to define onfinished handler';
    div.data('onfinished-callback', function() {docmill.compile(onfinished);});
    
    var skip=kwds.skip;
    if (!skip) skip = 500;
    var check, interval;
    check = function () {docmill.compile(onchecking);};
    interval = window.setInterval(check, skip);
    div.data('checking-interval', interval); div.data('skip', skip);

    // interior: status text and progressbar
    var status = tag('div'); status.addClass('status');
    div.append(status);
    if (kwds.status) status.text(kwds.status);

    var pbar = tag('div'); pbar.addClass('pbar');
    div.append(pbar);
    pbar.progressbar();
    if (kwds.percentage) pbar.progressbar('value', kwds.percentage);

    return div.lubanElement('progressbar');
  };
  //  object
  widgets.progressbar = function(elem) {
    this.super1 = widgets.base;
    this.super1(elem);
  };
  widgets.progressbar.prototype = new widgets.base;
  widgets.progressbar.prototype.setAttribute = function(attrs) {
    var je = this.jqueryelem;
    if (attrs.status) {
      var statusdiv = je.children('.status'); 
      statusdiv.text(attrs.status);
    }
    if (attrs.percentage) {
      var pbardiv = je.children('.pbar'); 
      pbardiv.progressbar('value', attrs.percentage);
      
      if (attrs.percentage>=100) {
	var interval = je.data('checking-interval');
	if (interval) {
	  clearInterval(interval); je.data('checking-interval', null);
	  var id = je.attr('id');
	  var f = function () {
	    var callback = je.data('onfinished-callback')
	    callback();
	  };
	  var skip = je.data('skip');
	  setTimeout(f, skip+200);
	}
      }
    }
  };
  widgets.progressbar.prototype.cancel = function(attrs) {
    var je = this.jqueryelem;
    clearInterval(je.data('checking-interval')); je.data('checking-interval', null);
  };


  // tabs
  //  factory
  ef.tabs = function(kwds, docmill) {

    var div = tag('div', {id: kwds.id} );

    var Class = kwds.Class;
    if (Class) div.addClass(Class);

    var ul = tag('ul');
    div.append(ul);

    var onclick = kwds.onclick;
    if (onclick) 
      div.click( function() { docmill.compile(onclick); return false; } );

    return div.lubanElement('tabs');
  };
  //  object
  widgets.tabs = function(elem) {
    this.super1 = widgets.base;
    this.super1(elem);
  };
  widgets.tabs.prototype = new widgets.base;
  widgets.tabs.prototype.createTab = function(tab, docmill) {
    var id = tab.id;
    var url = '#'+id;
    var label = tab.label;
    this._je.tabs('add', url, label);
    var tabdiv = $(url);

    var Class = tab.Class;
    if (Class) tabdiv.addClass(Class);
    
    // click
    var onclick = tab.onclick;
    if (onclick) {
      var callback = function() { 
	docmill.compile(onclick); return false;
      };
      tabdiv.bind('tabpanelshow', callback);
      tabdiv.click(callback);
    }

    // 
    tabdiv.bind('tabpanelshow', function() {$(this).trigger('resize');});

    // the last li corresponds to this tab
    var lis = this._je.children('ul').children('li');
    var lastli = $(lis[lis.length-1]);
    lastli.attr('target', id);
    
    return tabdiv.lubanElement('tab');
  };
  widgets.tabs.prototype.add = function (subelem) {
    throw "should not reach here";
  };


  // tab
  //  tab does not have standalone factory. only created from tabs.
  //  object
  widgets.tab = function(elem) {
    this.super1 = widgets.base;
    this.super1(elem);
  };
  widgets.tab.prototype = new widgets.base;
  widgets.tab.prototype.destroy = function() {
    var div = this._je;
    var parent = div.parent();
    var divs = parent.children('div');
    index = divs.index(div);
    parent.tabs('remove', index);
  };
  widgets.tab.prototype.setAttribute = function (attrs) {
    var div = this._je;
    var parent = div.parent();
    var ul = parent.children('ul');
    var id = div.attr('id');
    
    var label = attrs.label;
    if (label) {
      var li = ul.find("li[target='"+id+"']");
      span = li.find('span');
      span.text(label);
    }
  };


  // accordion
  //  factory
  ef.accordion = function(kwds, docmill, opts) {
    var plain_element = false;
    if (opts != null) {
      plain_element = opts.plain_element;
    }
    
    var div = tag('div', {id: kwds.id} );
    
    if (!plain_element) 
      div.accordion(ef.accordion.defaultopts);

    var onchange = kwds.onchange;
    if (onchange != null && onchange != '') {
      div.data('onchange-func', function(evt, ui) {
	  var postfixlen = 'label'.length;

	  var oldsection = ui.oldHeader.attr('id');
	  oldsection = oldsection.substr(0,oldsection.length-postfixlen);

	  var newsection = ui.newHeader.attr('id');
	  newsection = newsection.substr(0,newsection.length-postfixlen);

	  div.data('changed-data', {
	      'oldsection': oldsection, 'newsection': newsection});

	  docmill.compile(onchange); 

	  return false;

	});
      if (!plain_element)
	div.bind('accordionchange', div.data('onchange-func'));
    }

    return div.lubanElement('accordion');
  };
  ef.accordion.defaultopts = {
    //fillSpace: true,
  autoHeight: false
  };
  //  object
  widgets.accordion = function(elem) {
    this.super1 = widgets.base;
    this.super1(elem);
  };
  widgets.accordion.prototype = new widgets.base;
  widgets.accordion.prototype.createSection = function(section, opts) {
    
    var plain_element = false;
    if (opts != null) {
      plain_element = opts.plain_element;
    }
    
    var id = section.id;
    var label = section.label;
    
    var h3 = tag('h3', {'id': id+'label'});
    var a = tag('a'); h3.append(a);
    a.text(label);
    this._je.append(h3);
    
    var div = tag('div', {'id': id});
    this._je.append(div);
    
    // recreate the accordion
    if (!plain_element) {
      this._je.accordion('destroy');
      this._je.accordion(ef.accordion.defaultopts);

      // select the last section
      this._je.accordion('activate', $(h3));

    var onchange = div.data('onchange-func');
    if (onchange != null)
      this._je.bind('accordionchange', onchange);
    }

    return div.lubanElement('accordionsection');
  };
  widgets.accordion.prototype.removeSection = function(id) {
    this._je.accordion('destroy');
    var h3 = $('#'+id+'label');
    var div = $('#'+id);
    h3.remove();
    div.remove();
    
    // recreate the accordion
    this._je.accordion(ef.accordion.defaultopts);

    // all headers
    var headers = this._je.children('.ui-accordion-header');
    if (headers.length>1) {
      // select the first section
      this._je.accordion('activate', headers[0]);
    }
    
    //
    var onchange = div.data('onchange-func');
    if (onchange != null)
      this._je.bind('accordionchange', onchange);
    
    return;
  };
  widgets.accordion.prototype.add = function (subelem) {
    throw "should not reach here";
  };


  // accordionsection
  //  accordionsection does not have standalone factory. only created from accordion.
  //  object
  widgets.accordionsection = function(elem) {
    this.super1 = widgets.base;
    this.super1(elem);
  };
  widgets.accordionsection.prototype = new widgets.base;
  widgets.accordionsection.prototype.setAttribute = function (attrs) {
    var label = attrs.label;
    if (label != null) {
      var id = this._je.attr('id');
      var labelid = id+'label';
      var labelje = $('#'+labelid);
      labelje.find('a').text(label);
    }
  };
  widgets.accordionsection.prototype.destroy = function() {
    var div = this._je;
    var id = div.attr('id');
    
    var labelh3 = $(id+'label');
    labelh3.remove();
    
    var parent = div.parent()
    div.remove();
    
    parent.accordion('destroy');
    parent.accordion(ef.accordion.defaultopts);
    
    var onchange = div.data('onchange-func');
    if (onchange != null)
      parent.bind('accordionchange', onchange);
  };


  // appmenubar
  //  factory
  ef.appmenubar = function(kwds, docmill) {
    var ul = tag('ul', {id: kwds.id});
    return ul.lubanElement('appmenubar');
  };
  //  object
  widgets.appmenubar = function(elem) {
    this.super1 = widgets.base;
    this.super1(elem);
  };
  widgets.appmenubar.prototype = new widgets.base;

  // appmenu
  //  factory
  ef.appmenu = function(kwds, docmill) {
    var li = tag('li', {id: kwds.id});
    var span = tag('span'); li.append(span);
    span.addClass('luban-appmenu');

    var icon = kwds.icon;
    if (icon != null && icon != '') {
      var img = tag('img', {src: luban.iconpath(icon)});
      span.append(img);
    }
    span.append(kwds.label);
    
    var ul = tag('ul', {id: kwds.id+'interior-container'});
    li.append(ul);
    
    return li.lubanElement('appmenu');
  };
  //  object
  widgets.appmenu = function(elem) {
    this.super1 = widgets.base;
    this.super1(elem);
  };
  widgets.appmenu.prototype = new widgets.base;
  widgets.appmenu.prototype.add = function(subelem) {
    var icid = this._je.attr('id')+'interior-container';
    var ul = $('#'+icid);
    ul.append(subelem._je);
  };


  // appmenuitem
  //  factory
  ef.appmenuitem = function(kwds, docmill) {
    var li = tag('li', {id: kwds.id});
    var span = tag('span'); li.append(span);
    span.addClass('luban-appmenuitem');

    var icon = kwds.icon;
    if (icon != null && icon != '') {
      var img = tag('img', {src: luban.iconpath(icon)});
      span.append(img);
    }
    span.append(kwds.label);

    var onclick = kwds.onclick;
    if (onclick != null && onclick != '')
      li.click(function() {docmill.compile(kwds.onclick); return false;});

    return li.lubanElement('appmenuitem');
  };
  //  object
  widgets.appmenuitem = function(elem) {
    this.super1 = widgets.base;
    this.super1(elem);
  };
  widgets.appmenuitem.prototype = new widgets.base;


  // toolbar
  //  factory
  ef.toolbar = function(kwds, docmill) {
    var div = tag('div', {id: kwds.id} ); div.addClass('luban-toolbar');
    var table = tag('table'); div.append(table);
    var tr = tag('tr'); table.append(tr);
    tr.addClass('luban-toolbar-interior-container');
    
    div.data('luban-element-type', 'toolbar');
    return div.lubanElement();
  };
  //  object
  widgets.toolbar = function(elem) {
    this.super1 = widgets.base;
    this.super1(elem);
  };
  widgets.toolbar.prototype = new widgets.base;
  widgets.toolbar.prototype.add = function (subelem) {
    var td = tag('td');
    td.append(subelem.jqueryelem);
    
    var div = this.jqueryelem;
    var tr = div.find('.luban-toolbar-interior-container');
    tr.append(td);
  };


  // toolbarspacer
  //  factory
  ef.toolbarspacer = function(kwds, docmill) {
    var div = tag('div', {id: kwds.id});
    div.addClass('luban-toolbarspacer');
    div.data('luban-element-type', 'toolbarspacer');
    return div.lubanElement();
  };
  //  object
  widgets.toolbarspacer = function(elem) {
    this.super1 = widgets.base;
    this.super1(elem);
  };
  widgets.toolbarspacer.prototype = new widgets.base;


  // portlet
  //  factory
  ef.portlet = function(kwds, docmill) {
    var vpaddiv = tag('div', {id: kwds.id});
    vpaddiv.addClass('visualPadding');
    var div = tag('div'); div.addClass('portlet');
    vpaddiv.append(div);
    
    title = kwds.title
    if (title != null && title != '') {
      h5 = tag('h5'); div.append(h5);
      div.append(h5); h5.text(title);
    }
    
    bodydiv = tag('div'); div.append(bodydiv);
    bodydiv.addClass('portletBody');
    
    vpaddiv.data('luban-element-type', 'portlet');
    return vpaddiv.lubanElement();
  };
  //  object
  widgets.portlet = function(elem) {
    this.super1 = widgets.base;
    this.super1(elem);
  };
  widgets.portlet.prototype = new widgets.base;
  widgets.portlet.prototype.add = function (subelem) {
    var container = this.jqueryelem.children('div.portlet');
    var body = container.children('div.portletBody');
    body.append(subelem.jqueryelem);
  };
  widgets.portlet.prototype.setSelectedItem = function (item) {
    var container = this.jqueryelem.children('div.portlet');
    var body = container.children('div.portletBody');
    body.children('div.portletitem-container').removeClass('selected');
    item.jqueryelem.addClass('selected');
  };

  ef.portletitem = function(kwds, docmill) {

    var visualpadding = tag('div', {id: kwds.id}); 
    visualpadding.addClass('portletitem-container');

    var containerdiv = tag('div'); visualpadding.append(containerdiv);
    containerdiv.addClass('portletContent');

    var a = tag('a', {'title':kwds.tip})
    a.addClass(kwds.Class);
    containerdiv.append(a);

    // add icon if exists
    icon = kwds.icon
    if (icon != null && icon != '') {
      img = tag('img', {height: 16, width: 16, src: luban.iconpath(icon)});
      img.addClass(kwds.Class+'Icon');
      a.append(img);
    }
     
    // text
    span = tag('span');
    span.addClass(kwds.Class+"Text");
    span.text(kwds.label);
    a.append(span);

    // callbacks
    var onclick = kwds.onclick; var id = kwds.id;
    if (onclick != null && onclick != '')
      a.click( function() { 
	  var item = $('#'+id).lubanElement();
	  var portlet = item.getParent();
	  portlet.setSelectedItem(item);
	  docmill.compile(kwds.onclick); return false; 
	} );

    // selected?
    if (kwds.selected) visualpadding.addClass('selected');
    return visualpadding.lubanElement('portletitem');
  };
  //  object
  widgets.portletitem = function(elem) {
    this.super1 = widgets.base;
    this.super1(elem);
  };
  widgets.portletitem.prototype = new widgets.base;
  widgets.portletitem.prototype.getParent = function() {
    // see portlet.add and portlet factory
    var item = this._je;
    var body = item.parent();
    var portlet = body.parent();
    var vpad = portlet.parent();
    return vpad.lubanElement();
  };


  // form
  ef.form = function (kwds, docmill) {
    var Class = kwds.Class;
    var id = kwds.id;
    var ret = tag('form', {'id': id});
    ret.addClass(Class);

    var fieldset = tag('fieldset'); ret.append(fieldset);
    
    var title = kwds.title;
    var legend = tag('legend'); fieldset.append(legend);
    if (title) 
      legend.text(title);
    else 
      legend.hide();

    var onsubmit = kwds.onsubmit;
    if (onsubmit != null && onsubmit != '') {
      ret.submit(function () { docmill.compile(onsubmit); return false; });
    }
    
    var onclick = kwds.onclick;
    if (onclick) 
      ret.click( function() { docmill.compile(onclick); return false; } );

    return ret.lubanElement('form');
  };
  widgets.form = function(elem) {
    this.super1 = widgets.base;
    this.super1(elem);
  };
  widgets.form.prototype = new widgets.base;  
  widgets.form.prototype.add = function (subelem) {
    this._je.find('fieldset').append(subelem._je);
  };
  widgets.form.prototype.clearErrors = function (subelem) {
    // clear error boxes in the form
    this._je.find('.error').hide().text('');
  };
  widgets.form.prototype.setAttribute = function(attrs) {
    var form = this._je;
    var fieldset = form.children('fieldset')

    var title = attrs.title;
    var legend = fieldset.children('legend');
    if (title) {
      legend.text(title).show();
    } else {
      if (title=='') 
	legend.text(title).hide();
    }

    var Class = attrs.Class;
    if (Class) {
      form.removeClass();
      form.addClass(Class);
    }
  };


  // formtextfield
  ef.formtextfield = function(kwds, docmill) {
    var div = formfield(kwds, docmill);
    var field = kwds;
    var args =  {
      'name': prependActor(field.name),
      'type': 'text',
      'value': field.value,
      'id': field.id+'+input'
    };

    var input = tag('input', args); div.append(input);
    
    var onchange = kwds.onchange;
    if (onchange) {
      input.change( function() { docmill.compile(onchange); return false; } );
    }

    return div.lubanElement('formtextfield');
  };
  widgets.formtextfield = function(elem) {
    this.super1 = widgets.base;
    this.super1(elem);
  };
  widgets.formtextfield.prototype = new widgets.base;  
  widgets.formtextfield.prototype.setAttribute = function (attrs) {
    var je = this._je;
    formfield_setAttribute(je, attrs);
    
    var value = attrs.value;
    if (value != null) {
      input = je.children('input');
      input.val(value);
    }
  };

  
  // formpasswordfield
  ef.formpasswordfield = function(kwds, docmill) {
    var div = formfield(kwds, docmill);
    var field = kwds;
    var args =  {
      'name': prependActor(field.name),
      'type': 'password',
      'value': field.value,
      'id': field.id+'+input'
    };

    var input = tag('input', args); div.append(input);
    
    var onchange = kwds.onchange;
    if (onchange) {
      input.change( function() { docmill.compile(onchange); return false; } );
    }

    return div.lubanElement('formpasswordfield');
  };
  widgets.formpasswordfield = function(elem) {
    this.super1 = widgets.base;
    this.super1(elem);
  };
  widgets.formpasswordfield.prototype = new widgets.base;  
  widgets.formpasswordfield.prototype.setAttribute = function (attrs) {
    var je = this._je;
    formfield_setAttribute(je, attrs);
    
    var value = attrs.value;
    if (value != null) {
      input = je.children('input');
      input.val(value);
    }
  };


  // formsubmitbutton
  ef.formsubmitbutton = function(kwds, docmill) {
    var div = tag('div', {'id': kwds.id});
    div.addClass('luban-submitbutton-container');

    kls = kwds.Class;
    if (kls) div.addClass(kls);

    var field = kwds;
    var args =  {
      'type': 'submit',
      'value': field.label,
    };

    var input = tag('input', args);
    div.append(input);
    
    var onchange = kwds.onchange;
    if (onchange) {
      input.change( function() { docmill.compile(onchange); return false; } );
    }

    return div.lubanElement('formsubmitbutton');
  };
  widgets.formsubmitbutton = function(elem) {
    this.super1 = widgets.base;
    this.super1(elem);
  };
  widgets.formsubmitbutton.prototype = new widgets.base;  
  widgets.formsubmitbutton.prototype.setAttribute = function (attrs) {
    var je = this._je;
    
    var label = attrs.label;
    if (label) {
      var input = je.children('input');
      input.attr('value', label);
    }
  };

  // formtextarea
  ef.formtextarea = function(kwds, docmill) {
    var div = formfield(kwds, docmill);
    var field = kwds;
    var args =  {
      'name': prependActor(field.name),
      'id': field.id+'+input'
    };
    if (field.readonly) {
      args['readonly'] = 'readonly';
    }

    var input = tag('textarea', args);  div.append(input);
    input.text(field.value);
    
    var onchange = kwds.onchange;
    if (onchange) {
      input.change( function() { docmill.compile(onchange); return false; } );
    }

    return div.lubanElement('formtextarea');
  };
  widgets.formtextarea = function(elem) {
    this.super1 = widgets.base;
    this.super1(elem);
  };
  widgets.formtextarea.prototype = new widgets.base;  
  widgets.formtextarea.prototype.setAttribute = function (attrs) {
    var je = this._je;
    formfield_setAttribute(je, attrs);
    
    var input = je.children('textarea');
    
    var readonly = attrs.readonly;
    if (readonly != null) {
      if (readonly) 
	input.attr('readonly', 'readonly');
      else
	input.removeAttr('readonly');
    }

    var value = attrs.value;
    if (value != null) {
      input.text(value);
    }

  };


  // formselectorfield
  ef.formselectorfield = function(kwds, docmill) {
    var div = formfield(kwds, docmill);
    var field = kwds;
    var args =  {
      'name': prependActor(field.name),
      'id': field.id+'+input'
    };

    var input = tag('select', args);  div.append(input);
    
    var onchange = kwds.onchange;
    if (onchange) {
      input.change( function() { docmill.compile(onchange); return false; } );
    }

    var selection = field.selection;
    for (var i in field.entries) {
      var entry = field.entries[i];
      var value = entry.value, description;
      if (value == null) {
	value = entry[0]; description = entry[1];
      } else {
	description = entry.description;
      }
      var args = {'value': value};
      if (value == selection) {
	args.selected = 1;
      }
      var o = tag('option', args); input.append(o);
      o.text(description);
    }
    
    return div.lubanElement('formselectorfield');
  };
  widgets.formselectorfield = function(elem) {
    this.super1 = widgets.base;
    this.super1(elem);
  };
  widgets.formselectorfield.prototype = new widgets.base;  
  widgets.formselectorfield.prototype.setAttribute = function (attrs) {
    var je = this._je;
    formfield_setAttribute(je, attrs);
  };
  widgets.formselectorfield.prototype.getSelection = function (attrs) {
    var je = this._je;
    return je.find(':selected').val();
  };


  // formradiobox
  ef.formradiobox = function(kwds, docmill) {
    var div = formfield(kwds, docmill);
    var field = kwds;
    var args =  {
      'name': prependActor(field.name),
      'type': 'radio',
      //'value': field.value,
      'id': field.id+'+input'
    };

    for (var i in field.entries) {
      var entry = field.entries[i];
      var value = entry.value, description;
      if (value == null) {
	value = entry[0]; description = entry[1];
      } else {
	description = entry.description;
      }   
      args['value']=value; 
      var input = tag('input', args);
      var text = tag('label'); text.text(description)
      div.append(input);
      div.append(text);
    }

    // XXX: onchange?
    return div.lubanElement('formradiobox');
  };

  widgets.formradiobox = function(elem) {
    this.super1 = widgets.base;
    this.super1(elem);
  };
  widgets.formradiobox.prototype = new widgets.base;  
  widgets.formradiobox.prototype.setAttribute = function (attrs) {
    var je = this._je;
    formfield_setAttribute(je, attrs);
    
    var value = attrs.value;
    if (value != null) {
      input = je.children('input');
      input.val(value);
    }
  };


  // formfield
  function formfield( kwds, docmill ) {
    var id = kwds.id;
    
    var div = tag('div', {'id': id});
    div.addClass('formfield');
    if (kwds.Class) div.addClass(kwds.Class);

    var div1 = tag('div'); div.append(div1); div1.addClass('label-container');

    if (kwds.required) {
      var span = tag('span'); span.addClass("formfieldRequired"); div1.append(span);
      span.text('&nbsp;');
    }

    var labeldiv = tag('div'); div1.append(labeldiv);
    var label = tag('label', {'for': id}); labeldiv.append(label);
    if (kwds.label) {
      label.text(kwds.label);
    } else {
      labeldiv.hide();
    }

    
    var error = kwds.error;
    if (error == null) error = '';
    var errorid = '%s-error' % id;
    var errordiv = tag('div', {id: errorid}); div.append(errordiv);
    errordiv.addClass('error');
    if (error)
      errordiv.text(error);
    else
      errordiv.hide();

    var help = kwds.help;
    if (help == null) help = '';
    helpdiv = tag('div'); 
    helpdiv.addClass('formfieldHelp'); helpdiv.addClass('help');
    div.append(helpdiv);
    if (help)
      helpdiv.text(help);
    else
      helpdiv.hide();
    
    var onclick = kwds.onclick;
    if (onclick) 
      div.click( function() { docmill.compile(onclick); return false; } );

    return div;
  }

  function formfield_setAttribute(formfield, attrs) {
    var Class = attrs.Class;
    if (Class) {
      formfield.removeClass();
      formfield.addClass('formfield');
      formfield.addClass(Class);
    }

    var helpdiv = formfield.children('div.help');
    var help = attrs.help;
    if (help != null) {
      if (help) 
	helpdiv.text(help).show();
      else
	helpdiv.hide();
    }

    var errordiv = formfield.children('div.error');
    var error = attrs.error;
    if (error != null) {
      if (error) 
	errordiv.text(error).show();
      else
	errordiv.hide();
    }

    var labelcontainerdiv = formfield.children('div.label-container');
    var labeldiv = labelcontainerdiv.children('div');
    var label = attrs.label;
    if (label != null) {
      if (label) 
	labeldiv.text(label).show();
      else
	labeldiv.hide();
    }
  }

  
  // treeview
  //  factory
  ef.treeview = function(kwds, docmill) {
    var div = tag('div', {id: kwds.id});

//     var labeldiv = tag('div'); labeldiv.addClass('treeview-label');
//     div.append(labeldiv);
//     if (kwds.label) labeldiv.text(kwds.label);

//     var onclick = kwds.onclick;
//     if (onclick) {
//       div.click(function() {
// 	  luban.docmill.compile(onclick); 
// 	  return false;
// 	});
//     }
    
    var ul = tag('ul');  div.append(ul);
    ul.addClass('treeview-container');

    var onnodemoving = kwds.onnodemoving;
    if (onnodemoving) {
      var f = function(node, refnode, type, tree) {
	var data =  {
	  "node": $(node).attr('id'),
	  "refnode": $(refnode).attr('id'),
	  "position-type": type
	};
	div.data('nodemoving-data', data);
	
	docmill.compile(onnodemoving); 

	return false;
      };
      div.data('onnodemoving-func', f);
    }

    //var ul = tag('ul', {id: kwds.id});
    //ul.addClass('filetree');
    //return ul.lubanElement('treeview');
    return div.lubanElement('treeview');
  };
  //  object
  widgets.treeview = function(elem) {
    this.super1 = widgets.base;
    this.super1(elem);
  };
  widgets.treeview.prototype = new widgets.base;
  widgets.treeview.prototype.destroy = function() {
    var je = this._je;
    var id = je.attr('id');
    var ref = $.tree_reference(id);
    ref.destroy();
    je.remove();
  };
  widgets.treeview.prototype.add = function(subelem) {
    var div = this._je;
    var ul = div.find('.treeview-container');
    ul.append(subelem._je);
  };
  // set the root node
  // root: a luban TreeViewBranch instance
  widgets.treeview.prototype.setRoot = function(root) {
    var div = this._je;
    var ul = div.find('.treeview-container');
    var lis = ul.children('li');
    if (lis.length>0) 
      throw "there is already a root";
    ul.append(root._je);
  };
  widgets.treeview.prototype.removeNode = function(node) {
    var ref = this._jstreeRef();
    
    if (node) {
      ref.remove(node._je);
    } else {
      ref.remove();
    }

    ref.refresh();
  };
  widgets.treeview.prototype.getSelection = function(node) {
    var ref = this._jstreeRef();
    var selected = ref.selected;
    var ret = selected.lubanElement();
    return ret;
  };    
  widgets.treeview.prototype.closeAll = function() {
    var je = this._je;
    var id = je.attr('id');
    var ref = $.tree_reference(id);
    ref.close_all();
  };
  widgets.treeview.prototype.openBranch = function(branch) {
    var je = this._je;
    var id = je.attr('id');
    var ref = $.tree_reference(id);
    ref.open_branch(branch._je);
  };
  widgets.treeview.prototype.closeBranch = function(branch) {
    var je = this._je;
    var id = je.attr('id');
    var ref = $.tree_reference(id);
    ref.close_branch(branch._je);
  };
  widgets.treeview.prototype.selectNode = function(node) {
    var je = this._je;
    var id = je.attr('id');
    var ref = $.tree_reference(id);
    ref.select_branch(node._je);
  };
  widgets.treeview.prototype.addNode = function(referencenode, newnode, position) {
    var ref = this._jstreeRef();
    var obj = {
      'attributes': {id: newnode.id},
      'data': {title: newnode.label},
      'state': 'open',
      'children': []
    };
    var refnode = referencenode._je;
    ref.create(obj, refnode);
    ref.refresh();
    
    var li = $('#'+newnode.id);
    
    var onclick = newnode.onclick;
    if (onclick != null && onclick != '') {
      li.children('a').click(function() {
	  luban.docmill.compile(onclick); 
	  return true;
	});
    }

    return li.lubanElement(newnode.type);
  };
  // implementation
  widgets.treeview.prototype._jstreeRef = function() {
    var je = this._je;
    var id = je.attr('id');
    return $.tree_reference(id);
  };

  // treeviewbranch
  //  factory
  ef.treeviewbranch = function(kwds, docmill) {
    var li = tag('li', {id: kwds.id});

    var akwds = {};
    var icon = kwds.icon;
    if (icon != null && icon != '') {
      var iconpath = luban.iconpath(icon);
      akwds.style = 'background-image:url("'+iconpath+'");';
    }
    var a = tag("a", akwds); li.append(a);
    a.text(kwds.label);
    
    var ul = tag('ul'); ul.addClass('treeviewbranch-interior-container');
    li.append(ul);
    
    var onclick = kwds.onclick;
    if (onclick != null && onclick != '') 
      a.click(function() {
	  docmill.compile(onclick);
	  return true;
	});

    return li.lubanElement('treeviewbranch');
  };
  //  object
  widgets.treeviewbranch = function(elem) {
    this.super1 = widgets.base;
    this.super1(elem);
  };
  widgets.treeviewbranch.prototype = new widgets.base;
  widgets.treeviewbranch.prototype.add = function(subelem) {
    var ul = this._je.children('ul.treeviewbranch-interior-container');
    ul.append(subelem._je);
  };


  // treeviewleaf
  //  factory
  ef.treeviewleaf = function(kwds, docmill) {
    var li = tag('li', {id: kwds.id});

    var akwds = {};
    var icon = kwds.icon;
    if (icon != null && icon != '') {
      var iconpath = luban.iconpath(icon);
      akwds.style = 'background-image:url("'+iconpath+'");';
    }
    var a = tag("a", akwds); li.append(a);
    a.text(kwds.label);

    var onclick = kwds.onclick;
    if (onclick != null && onclick != '') 
      a.click(function() {docmill.compile(onclick); return true;});
    return li.lubanElement('treeviewleaf');
  };
  //  object
  widgets.treeviewleaf = function(elem) {
    this.super1 = widgets.base;
    this.super1(elem);
  };
  widgets.treeviewleaf.prototype = new widgets.base;




  // helpers
  function tag(name, kwds) {

    var assignments = [];
    
    for (var key in kwds) {
      var value = kwds[key];
      assignments.push( key + '=' + '"' + value + '"' );
    }
    
    var  s = "<" + name + ' ' + assignments.join(' ') + ">" + "</"+name+">";
    
    return $(s);
  };

  // prepend 'actor.' to keys
  function prependActor(s) {
    return 'actor.'+s;
  }

 })(luban, jQuery);


// End of file
