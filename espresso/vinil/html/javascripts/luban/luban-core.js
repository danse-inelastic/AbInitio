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

// namespace luban
luban = {
  'elementFactory': {},
  'widgets': {},
  'configuration': {
    'images_base': '/images',
    'icons_base': '/images/icons',
  }
};



(function(luban, $) {


  $.fn.lubanElement = function (type) {
    if (type==null)
      type = this.data('luban-element-type');
    else
      this.data('luban-element-type', type);
    var factory = luban.widgets[type];
    return new factory(this);
  };

  // widget base
  luban.widgets.base = function (elem) {
    this.jqueryelem = elem;
    this._je = this.jqueryelem;
  };
  luban.widgets.base.prototype = {
    'type': function () {
      return this.jqueryelem.data('luban-element-type');
    },
    'add': function (subelem) {
      if (typeof(subelem) == 'string') {
	this.jqueryelem.append(subelem);
      } else {
	this.jqueryelem.append(subelem.jqueryelem);
      }
    },
    'destroy': function() {
      this.jqueryelem.remove();
    },
    'setAttribute': function(args) {
      throw 'widgets.base.setAttribute:' + this.type() + ' notimplementederror';
    },
    
    // retrieve data related to the specified event
    'getEventData': function (event) {
      return this.jqueryelem.data(event+'-data');
    },

    // empty my content
    'empty': function (event) {
      this.jqueryelem.empty();
      //throw 'widgets.base.empty:' + this.type() +' notimplementederror';
    },

    // addClass
    'addClass': function(Class) {
      this.jqueryelem.addClass(Class);
    },

    // removeClass
    'removeClass': function(Class) {
      if (this.jqueryelem.hasClass(Class)) {
	this.jqueryelem.removeClass(Class);
      }
    },

  };

  luban.iconpath = function(filename) {
    return luban.configuration.icons_base+'/'+filename;
  };

  luban.imagepath = function(filename) {
    return luban.configuration.images_base+'/'+filename;
  };

 })(luban, jQuery);


// End of file
