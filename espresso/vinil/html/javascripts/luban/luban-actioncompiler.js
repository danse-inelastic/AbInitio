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
//    * luban-controller.js


(function(luban, $) {

  luban.actioncompiler = function (docmill) {
    if (docmill == null) docmill = new luban.documentmill;
    this.docmill = docmill;
  };
  
  luban.actioncompiler.prototype = {
    
    'compile': function(actions) {
      // check if it is an action of a list of actions
      var type = actions.type;
      if (type != null) return this.compile1(actions);
      
      // 
      for (var i in actions) {
	var action = actions[i];
	this.compile1(action);
      }
    },

    'compile1': function (action) {
      return this.dispatch(action);
    },

    'dispatch': function (action) {
      var type = action.type;
      code = 'this.on'+type+"(action)";
      return eval(code);
    },

    'onloading': function(action) {
      var params = {};
      for (var key in action.params) {
	var value = action.params[key];
	if (typeof(value) == 'object') {
	  value = this.dispatch(value);
	}
	params[key] = value;
      }
      var kwds = {
	'actor': action.actor,
	'routine': action.routine,
	'data': params
      };
      var C = luban.Controller;
      return C.load(kwds);
    },

    'onsubmission': function(action) {
      var form = action.form;
      form = this.dispatch(form);
      
      var kwds = {
	'actor': action.actor,
	'routine': action.routine,
	'data': action.params
      };

      form.clearErrors();
      var C = luban.Controller; C.submit(form._je, kwds);
    },

    'onselectbyid': function(action) {
      var id = action.id;
      if (!id) {
        // assume that no id means the page
        return $('div.luban-page').lubanElement();
      }
      var je = $('#'+id);
      if (je.length == 0) {
	throw "not such element: id="+id;
      }
      return je.lubanElement();
    },

    'onreplacecontent': function(action) {
      var e = action.element;
      var element = this.dispatch(e);
      element.empty();
      
      var newdoc = action.newcontent;
      this.docmill.render(newdoc, element);

      element.jqueryelem.trigger('resize');
    },

    'onremovecontent': function(action) {
      var e = action.element;
      var element = this.dispatch(e);
      element.empty();
    },

    'onsimpleaction': function(action) {
      var name = action.actionname;
      return eval('this.on'+name+'(action.params)');
    },

    'onsimpleelementaction': function(action) {
      var element = action.element;
      element = this.dispatch(element);
      switch(action.actionname) {

      case 'show':
	return element.show();

      case 'hide':
	return element.hide();

      case 'destroy':
	return element.destroy();

      case 'addClass':
	return element.addClass(action.params.Class);

      case 'removeClass':
	return element.removeClass(action.params.Class);

      case 'setAttribute':
	return element.setAttribute(action.params);

      default:
	var etype = element._je.data('luban-element-type');
	var method = 'on'+etype+action.actionname;
	return eval('this.'+method+'(action);');
      }
    },

    'onformfieldshowerrormessage': function(action) {
      return this.onformfieldshowerror(action);
    },
    'onformfieldshowerror': function(action) {
      var params = action.params;
      var message = params.message;
      element = this.dispatch(action.element);
      return element._je.find('.error').text(message).show('normal');
    },

    'onformtextfieldshowerrormessage': function(action) {
      this.onformfieldshowerrormessage(action);
    },
    'onformtextfieldshowerror': function(action) {
      this.onformfieldshowerror(action);
    },
    'onformpasswordfieldshowerrormessage': function(action) {
      this.onformfieldshowerrormessage(action);
    },
    'onformpasswordfieldshowerror': function(action) {
      this.onformfieldshowerror(action);
    },
    'onformselectorfieldshowerrormessage': function(action) {
      this.onformfieldshowerrormessage(action);
    },
    'onformselectorfieldshowerror': function(action) {
      this.onformfieldshowerror(action);
    },
    'onformtextareashowerrormessage': function(action) {
      this.onformfieldshowerrormessage(action);
    },
    'onformtextareashowerror': function(action) {
      this.onformfieldshowerror(action);
    },
    'onformselectorfieldgetSelection': function(action) {
      var element = this.dispatch(action.element);
      return element.getSelection();
    },
    'onformtextfieldgetValue': function(action) {
      var element = this.dispatch(action.element);
      return element._je.find('input').val();
    },

    'onprogressbarcancel': function(action) {
      var element = this.dispatch(action.element);
      element.cancel();
    },

    'oncredentialremoval': function(action) {
      luban.Controller.credential = null;
    },

    'onappendelement': function(action) {
      var container = this.compile(action.container);
      this.docmill.render(action.element, container);
    },

    'onnotification': function(action) {
      var element = this.dispatch(action.element);
      var event = action.event;
      var kwds = {
	'actor': action.actor,
	'routine': action.routine,
	'data': action.params
      };
      var C = luban.Controller;
      return C.notify(element, event, kwds);
    },

    'ontreeviewsetroot': function(action) {
      var treeview = action.treeview;
      treeview = this.dispatch(treeview);
      
      var root = this.docmill.render(action.root);
      treeview.setRoot(root);
    },

    'ontreeviewaddbranch': function(action) {
      var treeview = action.treeview;
      treeview = this.dispatch(treeview);
      
      var referencenode = this.dispatch(action.referencenode);
      
      treeview.addNode(referencenode, action.newnode, action.position);
    },

    'ontreeviewcloseAll': function(action) {
      var treeview = this.dispatch(action.element);

      treeview.closeAll();
    },


    'ontreeviewopen': function(action) {
      var treeview = this.dispatch(action.element);

      treeview.open(action.params.branch);
    },

    'ontreeviewclose': function(action) {
      var treeview = this.dispatch(action.element);

      treeview.close(action.params.branch);
    },

    'ontreeviewselectnode': function(action) {
      var treeview = this.dispatch(action.treeview);
      var node = this.dispatch(action.node);

      treeview.selectNode(node);
    },

    'ontreeviewremovenode': function(action) {
      var treeview = action.treeview;
      treeview = this.dispatch(treeview);
      
      var node = this.dispatch(action.node);
      treeview.removeNode(node);
    },


    'ontreeviewgetSelection': function(action) {
      var treeview = this.dispatch(action.element);
      return treeview.getSelection();
    },


    'onaccordioncreateSection': function(action) {
      var accordion = this.dispatch(action.element);
      accordion.createSection(action.params);
    },

    'onaccordionremoveSection': function(action) {
      var accordion = this.dispatch(action.element);
      accordion.removeSection(action.params.id);
    },

    'onalert': function(params) {
      alert(params.message);
    }

  };

 })(luban, jQuery);


// End of file
