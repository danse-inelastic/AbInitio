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
      var kwds = {
	'actor': action.actor,
	'routine': action.routine,
	'data': action.params,
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
	'data': action.params,
      };

      var C = luban.Controller;
      C.clearFormErrorAlerts(form);
      C.submit(form, kwds);
    },

    'onselectbyid': function(action) {
      var id = action.id;
      return $('#'+id);
    },

    'onreplacecontent': function(action) {
      var e = action.element;
      var element = this.dispatch(e);
      element.lubanElement().empty();
      
      var newdoc = action.newcontent;
      this.docmill.render(newdoc, element.lubanElement());
    },

    'onremovecontent': function(action) {
      var e = action.element;
      var element = this.dispatch(e);
      element.lubanElement().empty();
    },

    'onsimpleaction': function(action) {
      var name = action.actionname;
      return eval('this.on'+name+'(action.params)');
    },

    'onsimpleelementaction': function(action) {
      var element = action.element;
      element = this.dispatch(element);
      switch(action.actionname) {

      case 'destroy':
	return element.lubanElement().destroy();

      case 'addClass':
	return element.lubanElement().addClass(action.params.Class);

      case 'removeClass':
	return element.lubanElement().removeClass(action.params.Class);

      case 'setAttribute':
	return element.lubanElement().setAttribute(action.params);

      default:
	var etype = element.data('luban-element-type');
	if (etype==null) throw 'no element type';
	
	var method = 'on'+etype+action.actionname;
	eval('this.'+method+'(action);');
      }
    },

    'onformfieldshowerrormessage': function(action) {
      var params = action.params;
      var message = params.message;
      element = this.dispatch(action.element);
      return element.find('.error').text(message).show('normal');
    },

    'onformtextfieldshowerrormessage': function(action) {
      this.onformfieldshowerrormessage(action);
    },
    'onformpasswordfieldshowerrormessage': function(action) {
      this.onformfieldshowerrormessage(action);
    },
    'onformselectorfieldshowerrormessage': function(action) {
      this.onformfieldshowerrormessage(action);
    },
    'onformtextareashowerrormessage': function(action) {
      this.onformfieldshowerrormessage(action);
    },


    'oncredentialremoval': function(action) {
      luban.Controller.credential = null;
    },

    
    'onappendelement': function(action) {
      var container = this.compile(action.container);
      this.docmill.render(action.element, container.lubanElement());
    },


    'onnotification': function(action) {
      var element = this.dispatch(action.element);
      element = element.lubanElement();
      var event = action.event;
      var kwds = {
	'actor': action.actor,
	'routine': action.routine,
	'data': action.params,
      };
      var C = luban.Controller;
      return C.notify(element, event, kwds);
    },

    'ontreeviewsetroot': function(action) {
      var treeview = action.treeview;
      treeview = this.dispatch(treeview);
      treeview = treeview.lubanElement();
      
      var root = this.docmill.render(action.root);
      treeview.setRoot(root);
    },

    'ontreeviewaddbranch': function(action) {
      var treeview = action.treeview;
      treeview = this.dispatch(treeview);
      treeview = treeview.lubanElement();
      
      treeview.addNode(action.referencenode, action.newnode, action.position);
    },

    'ontreeviewremovenode': function(action) {
      var treeview = action.treeview;
      treeview = this.dispatch(treeview);
      treeview = treeview.lubanElement();
      
      var node = action.node;
      node = this.dispatch(node);
      node = node.lubanElement();

      treeview.removeNode(node);
    },


    'ontreeviewcloseAll': function(action) {
      var treeview = this.dispatch(action.element);
      treeview = treeview.lubanElement();

      treeview.closeAll();
    },


    'ontreeviewopen': function(action) {
      var treeview = this.dispatch(action.element);
      treeview = treeview.lubanElement();

      treeview.open(action.params.branch);
    },

    'ontreeviewclose': function(action) {
      var treeview = this.dispatch(action.element);
      treeview = treeview.lubanElement();

      treeview.close(action.params.branch);
    },

    'ontreeviewselect': function(action) {
      var treeview = this.dispatch(action.element);
      treeview = treeview.lubanElement();

      treeview.select(action.params.branch);
    },

    'onaccordioncreateSection': function(action) {
      var accordion = this.dispatch(action.element);
      accordion = accordion.lubanElement();
      accordion.createSection(action.params);
    },

    'onaccordionremoveSection': function(action) {
      var accordion = this.dispatch(action.element);
      accordion = accordion.lubanElement();
      accordion.removeSection(action.params.id);
    },

    'onalert': function(params) {
      alert(params.message);
    },

  };

 })(luban, jQuery);


// End of file
