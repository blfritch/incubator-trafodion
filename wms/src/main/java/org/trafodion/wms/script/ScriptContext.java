/**
 *(C) Copyright 2015 Hewlett-Packard Development Company, L.P.
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 *     http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 */
package org.trafodion.wms.script;

import java.util.*;
import java.io.IOException;

import org.apache.commons.logging.Log;
import org.apache.commons.logging.LogFactory;

import org.trafodion.wms.Constants;

public final class ScriptContext {
	private static final Log LOG = LogFactory.getLog(ScriptContext.class);
	private String scriptName;
	private String hostName;
	private String command;
	private boolean debug=false;
	private int exitCode=0;
	private StringBuilder stdOut=new StringBuilder();
	private StringBuilder stdErr=new StringBuilder();

	public ScriptContext(){
	}
	public void setDebug(boolean value){
		debug = value;
	}
	public boolean isDebug(){
		return debug;
	}
	public void setScriptName(String value){
		this.scriptName=value;
	}
	public String getScriptName(){
		return scriptName;
	}
	public void setHostName(String value){
		this.hostName=value;
	}
	public String getHostName(){
		return hostName;
	}
	public void setCommand(String value){
		this.command=value;
	}
	public String getCommand(){
		return command;
	}
	public void setExitCode(int value){
		this.exitCode=value;
	}
	public int getExitCode(){
		return exitCode;
	}
	public StringBuilder getStdOut(){
		return stdOut;
	}
	public StringBuilder getStdErr(){
		return stdErr;
	}
	public String toString(){
		StringBuilder sb = new StringBuilder();
		sb.append("\nScript name[" + getScriptName() + "]");
		sb.append("\nCommand[" + getCommand() + "]");
		sb.append("\nHost name[" + getHostName() + "]");
		sb.append("\nExit code[" + getExitCode() + "]");
		sb.append("\nStandard out[" + getStdOut().toString() + "]");
		sb.append("\nStandard error[" + getStdErr().toString() + "]");
		return sb.toString();
	}
}
